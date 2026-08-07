"""Tests for pyRiverBed.

Uses the standard library's ``unittest`` rather than pytest, so that running
the suite needs nothing beyond the package's own dependencies::

    python -m unittest discover -s tests -v
"""

import math
import re
import tempfile
import unittest
from pathlib import Path

import numpy as np

from pyriverbed import bed, geometry, mesh, migration, planform
from pyriverbed.config import (ChannelConfig, ChuteCutoffConfig, Config,
                               ConfigError, CurvatureConfig, FlipConfig,
                               KinoshitaConfig, MigrationConfig,
                               NeckCutoffConfig, default_config, dump_ini,
                               read_ini, read_legacy_steering)
from pyriverbed.model import RiverBedModel
from pyriverbed.planform import Centerline


def sine_centerline(n=2000, spacing=1.0, wavelength=400.0, amplitude=60.0):
    """A simple sinusoidal centerline for the cutoff tests."""
    s = np.arange(n) * spacing
    x = s
    y = amplitude * np.sin(2 * np.pi * s / wavelength)
    line_s = planform.arc_length(x, y)
    cur, theta = planform.curvature(line_s, x, y)
    return Centerline(x=x, y=y, s=line_s, curvature=cur, theta=theta)


class TestConfig(unittest.TestCase):
    def test_defaults_validate(self):
        self.assertIsInstance(default_config(), Config)

    def test_ini_roundtrip(self):
        config = default_config()
        config.channel.width = 12.5
        config.chute_cutoff.enabled = True
        config.chute_cutoff.entrance = "inflection"
        config.migration.seed = 99
        self.assertEqual(read_ini(dump_ini(config)).to_dict(),
                         config.to_dict())

    def test_seed_none_roundtrip(self):
        config = default_config()
        config.migration.seed = None
        self.assertIsNone(read_ini(dump_ini(config)).migration.seed)

    def test_booleans_accept_words(self):
        text = dump_ini(default_config()).replace("despike         = yes",
                                                  "despike         = true")
        self.assertTrue(read_ini(text).curvature.despike)

    def test_unknown_key_is_rejected(self):
        text = dump_ini(default_config()) + "\nnot_a_key = 3\n"
        with self.assertRaises(ConfigError):
            read_ini(text)

    def test_out_of_range_is_rejected(self):
        for section, key, value in (
            ("channel", "width", -1.0),
            ("channel", "transverse_slope_corrector", 2.0),
            ("chute_cutoff", "frequency", 1.5),
            ("chute_cutoff", "max_valley_angle", 120.0),
            ("chute_cutoff", "min_sinuosity", 0.5),
            ("curvature", "method", "bogus"),
            ("chute_cutoff", "entrance", "bogus"),
        ):
            with self.subTest(key=f"{section}.{key}"):
                config = default_config()
                setattr(getattr(config, section), key, value)
                with self.assertRaises(ConfigError):
                    config.validate()

    def test_derived_channel_quantities(self):
        channel = ChannelConfig(width=0.6, depth=0.15, n_offsets=10)
        self.assertAlmostEqual(channel.half_width, 0.3)
        self.assertAlmostEqual(channel.interval, 0.03)
        self.assertAlmostEqual(channel.beta, 2.0)
        beta = 2.0
        self.assertAlmostEqual(
            channel.scour_factor,
            3.8 * (1 + beta / 6.96 * math.exp(-6.96 / beta)))

    def test_smoothing_level_escalates_above_38(self):
        self.assertEqual(CurvatureConfig(smoothing_level=20).n_passes, 20)
        self.assertEqual(CurvatureConfig(smoothing_level=50).n_passes,
                         int(round(1.1 ** 50)))

    def test_legacy_steering_v11(self):
        values = ["jurua.txt", 2, 3, 10.0, 110.0, 0.03125, 0.00520833, 160.0,
                  8.0, 0.0, 0.03, 10, 1, 6.0, 1, 1, 1, 0, 0, 1, 3.0, 0.0, 0.01,
                  0.1, 86400.0, 1e-07, 50, 5000, 100, 24, 20, 1.0, 1, 0.2,
                  1000, 2, 2, 45, 8, 3]
        config = read_legacy_steering("\n".join(str(v) for v in values))
        self.assertEqual(config.mode, "centerline")
        self.assertEqual(config.centerline_file, "jurua.txt")
        self.assertAlmostEqual(config.channel.width, 160.0)
        self.assertAlmostEqual(config.lag.strength, 6.0)
        self.assertTrue(config.migration.enabled)
        self.assertTrue(config.chute_cutoff.enabled)
        self.assertEqual(config.chute_cutoff.entrance, "inflection")
        self.assertAlmostEqual(config.chute_cutoff.frequency, 0.2)
        self.assertEqual(config.chute_cutoff.span, 2)

    def test_legacy_steering_v10_has_chute_off(self):
        values = ["jurua.txt", 2, 3, 10.0, 110.0, 0.03125, 0.00520833, 160.0,
                  8.0, 0.0, 0.03, 10, 1, 6.0, 1, 1, 1, 0, 0, 1, 3.0, 0.0, 0.01,
                  0.1, 86400.0, 1e-07, 50, 5000, 100, 24, 20, 1.0]
        config = read_legacy_steering("\n".join(str(v) for v in values))
        self.assertFalse(config.chute_cutoff.enabled)

    def test_truncated_steering_is_rejected(self):
        with self.assertRaises(ConfigError):
            read_legacy_steering("f.txt\n1\n2\n3\n")


class TestPlanform(unittest.TestCase):
    def test_arc_length_of_straight_line(self):
        x = np.arange(5, dtype=float)
        s = planform.arc_length(x, np.zeros(5))
        np.testing.assert_allclose(s, x)

    def test_kinoshita_is_closed_and_periodic(self):
        line = planform.build_kinoshita(
            KinoshitaConfig(n_bends=3, arc_wavelength=10.0),
            ChannelConfig(ds=0.01))
        # Arc length must match the requested number of wavelengths.
        self.assertAlmostEqual(line.s[-2], 30.0, places=6)
        # A sine-generated curve returns to zero transverse offset each period.
        self.assertGreater(line.sinuosity, 1.0)

    def test_kinoshita_without_harmonics_is_sine_generated(self):
        line = planform.build_kinoshita(
            KinoshitaConfig(n_bends=2, arc_wavelength=10.0, skewness=0.0,
                            flatness=0.0, max_angular_amplitude=30.0),
            ChannelConfig(ds=0.005))
        # theta is an odd function about the half wavelength, so y returns to 0.
        n_period = int(10.0 / 0.005)
        self.assertAlmostEqual(line.y[n_period], 0.0, places=3)

    def test_resample_of_a_straight_line_is_exactly_uniform(self):
        x = np.array([0.0, 1.0, 4.0, 10.0])
        s, xn, yn = planform.resample(x, np.zeros(4), 0.25)
        spacing = np.diff(s)
        np.testing.assert_allclose(spacing, spacing[0], atol=1e-12)
        self.assertAlmostEqual(xn[0], 0.0)
        self.assertAlmostEqual(xn[-1], 10.0)

    def test_resample_is_near_uniform_and_keeps_ends(self):
        # Stations are equally spaced along the *input* arc length, so on a
        # curved polyline the resampled chords differ from each other by the
        # curvature-induced shortening only.
        s0 = np.linspace(0, 400, 500)
        x = s0
        y = 30 * np.sin(2 * np.pi * s0 / 200)
        s, xn, yn = planform.resample(x, y, 1.0)
        spacing = np.diff(s)
        np.testing.assert_allclose(spacing, np.mean(spacing), rtol=1e-3)
        self.assertAlmostEqual(xn[0], x[0])
        self.assertAlmostEqual(yn[-1], y[-1])

    def test_smooth_pins_the_ends(self):
        rng = np.random.default_rng(0)
        x = np.linspace(0, 10, 200)
        y = np.sin(x) + rng.normal(0, 0.05, x.size)
        _, xs, ys = planform.smooth(x, y, 10)
        self.assertAlmostEqual(xs[0], x[0])
        self.assertAlmostEqual(ys[-1], y[-1])
        # Smoothing must reduce the roughness.
        self.assertLess(np.std(np.diff(ys, 2)), np.std(np.diff(y, 2)))

    def test_curvature_of_a_circle(self):
        radius = 25.0
        angle = np.linspace(0, np.pi, 4000)
        x, y = radius * np.cos(angle), radius * np.sin(angle)
        s = planform.arc_length(x, y)
        for method in CurvatureConfig.METHODS:
            with self.subTest(method=method):
                cur, _ = planform.curvature(
                    s, x, y, CurvatureConfig(method=method, despike=False))
                interior = np.abs(cur[10:-10])
                np.testing.assert_allclose(interior, 1 / radius, rtol=2e-3)

    def test_curvature_sign_follows_turn_direction(self):
        angle = np.linspace(0, np.pi, 2000)
        # Increasing angle traverses the upper semicircle counterclockwise,
        # i.e. a left turn, which must give positive curvature.
        x, y = 20 * np.cos(angle), 20 * np.sin(angle)
        s = planform.arc_length(x, y)
        cur, _ = planform.curvature(s, x, y)
        self.assertGreater(cur[1000], 0)
        cur_mirrored, _ = planform.curvature(s, x, -y)
        self.assertLess(cur_mirrored[1000], 0)

    def test_curvature_survives_a_westward_reach(self):
        # A reach heading due west straddles the +/-pi branch cut of arctan2.
        # Without angle wrapping this produced an enormous false spike.
        angle = np.linspace(0.5 * np.pi, 1.5 * np.pi, 3000)
        x, y = 30 * np.cos(angle), 30 * np.sin(angle)
        s = planform.arc_length(x, y)
        cur, _ = planform.curvature(s, x, y, CurvatureConfig(despike=False))
        np.testing.assert_allclose(np.abs(cur[10:-10]), 1 / 30.0, rtol=2e-3)

    def test_despike_removes_a_single_spike(self):
        cur = np.full(21, 0.01)
        cur[10] = 5.0
        cleaned = planform.despike_curvature(cur)
        self.assertAlmostEqual(cleaned[10], 0.01)

    def test_phase_lag_preserves_mean_and_shifts_downstream(self):
        s = np.linspace(0, 100, 2000)
        cur = np.sin(2 * np.pi * s / 25)
        lagged = planform.phase_lag(cur, 120)
        self.assertAlmostEqual(np.mean(lagged[300:]), np.mean(cur[300:]),
                               places=2)
        # The lagged peak must sit downstream of the original peak.
        window = slice(400, 900)
        self.assertGreater(np.argmax(lagged[window]), np.argmax(cur[window]))

    def test_phase_lag_weights_sum_to_one(self):
        constant = np.ones(500)
        lagged = planform.phase_lag(constant, 60)
        np.testing.assert_allclose(lagged[100:], 1.0, atol=1e-12)

    def test_phase_lag_noop_for_tiny_window(self):
        cur = np.linspace(0, 1, 50)
        np.testing.assert_allclose(planform.phase_lag(cur, 1), cur)

    def test_extend_ends_adds_straight_tangent_reaches(self):
        x = np.linspace(0, 10, 101)
        y = np.zeros_like(x)
        xe, ye = planform.extend_ends(x, y, extension=5.0, spacing=0.1)
        self.assertAlmostEqual(xe[0], -5.0, places=6)
        self.assertAlmostEqual(xe[-1], 15.0, places=6)
        np.testing.assert_allclose(ye, 0.0, atol=1e-12)

    def test_load_centerline_rejects_one_column(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.txt"
            path.write_text("1\n2\n3\n")
            with self.assertRaises(ValueError):
                planform.load_centerline(path)


class TestBed(unittest.TestCase):
    def setUp(self):
        self.channel = ChannelConfig(width=10.0, depth=1.0, n_offsets=20,
                                     slope=0.0)
        self.flip = FlipConfig(transverse=False)

    def test_straight_channel_is_flat(self):
        s = np.linspace(0, 100, 200)
        result = bed.compute_bed(np.zeros_like(s), s, self.channel, self.flip)
        np.testing.assert_allclose(result.z, 0.0, atol=1e-6)
        np.testing.assert_allclose(result.centerline_depth, self.channel.depth,
                                   rtol=1e-6)

    def test_cross_sectional_area_is_conserved(self):
        # This is the constraint that fixes the centerline depth, so it is the
        # single most important invariant of the bed model.
        s = np.linspace(0, 500, 400)
        cur = 0.02 * np.sin(2 * np.pi * s / 120)
        result = bed.compute_bed(cur, s, self.channel, self.flip)
        depth_profile = self.channel.depth - result.z
        n = np.linspace(-self.channel.half_width, self.channel.half_width,
                        2 * self.channel.n_offsets + 1)
        area = np.trapezoid(depth_profile, n, axis=1)
        expected = self.channel.width * self.channel.depth
        np.testing.assert_allclose(area, expected, rtol=2e-3)

    def test_pool_sits_on_the_outer_bank(self):
        # Column 0 is the right bank and the last column is the left bank (see
        # geometry.build_point_cloud). A steady left turn has its outer bank on
        # the right, so the pool must be at column 0 and the bar at the last.
        s = np.linspace(0, 100, 50)
        cur = np.full_like(s, 0.02)            # a steady left turn
        section = bed.compute_bed(cur, s, self.channel, self.flip).z[25]
        self.assertLess(section[0] * section[-1], 0)
        self.assertEqual(int(np.argmin(section)), 0)
        self.assertEqual(int(np.argmax(section)), section.size - 1)

    def test_right_turn_puts_the_pool_on_the_left(self):
        s = np.linspace(0, 100, 50)
        cur = np.full_like(s, -0.02)           # a steady right turn
        section = bed.compute_bed(cur, s, self.channel, self.flip).z[25]
        self.assertEqual(int(np.argmin(section)), section.size - 1)
        self.assertEqual(int(np.argmax(section)), 0)

    def test_reversing_curvature_mirrors_the_bed(self):
        s = np.linspace(0, 100, 50)
        cur = np.full_like(s, 0.02)
        left = bed.compute_bed(cur, s, self.channel, self.flip)
        right = bed.compute_bed(-cur, s, self.channel, self.flip)
        np.testing.assert_allclose(left.z, right.z[:, ::-1], atol=1e-9)

    def test_transverse_slope_corrector_scales_relief(self):
        s = np.linspace(0, 100, 50)
        cur = np.full_like(s, 0.02)
        full = bed.compute_bed(cur, s, self.channel, self.flip).relief
        half = ChannelConfig(width=10.0, depth=1.0, n_offsets=20,
                             transverse_slope_corrector=0.5)
        self.assertLess(bed.compute_bed(cur, s, half, self.flip).relief, full)

    def test_flip_transverse_mirrors_the_bed(self):
        s = np.linspace(0, 100, 50)
        cur = np.full_like(s, 0.02)
        normal = bed.compute_bed(cur, s, self.channel, FlipConfig(
            transverse=False))
        flipped = bed.compute_bed(cur, s, self.channel, FlipConfig(
            transverse=True))
        np.testing.assert_allclose(normal.z, flipped.z[:, ::-1])

    def test_longitudinal_slope_raises_the_upstream_end(self):
        s = np.linspace(0, 1000, 100)
        channel = ChannelConfig(width=10.0, depth=1.0, n_offsets=5, slope=1e-3)
        result = bed.compute_bed(np.zeros_like(s), s, channel, self.flip)
        self.assertAlmostEqual(result.z[0, 2] - result.z[-1, 2], 1.0, places=3)


class TestGeometry(unittest.TestCase):
    def test_offset_of_a_straight_line(self):
        x = np.linspace(0, 10, 51)
        y = np.zeros_like(x)
        xo, yo = geometry.offset_polyline(x, y, 2.0)
        np.testing.assert_allclose(yo[:, 0], 2.0, atol=1e-9)
        np.testing.assert_allclose(yo[:, 1], -2.0, atol=1e-9)
        np.testing.assert_allclose(xo[:, 0], x, atol=1e-9)

    def test_offset_keeps_width_constant_around_a_bend(self):
        angle = np.linspace(0, np.pi / 2, 400)
        x, y = 50 * np.cos(angle), 50 * np.sin(angle)
        xo, yo = geometry.offset_polyline(x, y, 5.0)
        separation = np.hypot(xo[:, 0] - xo[:, 1], yo[:, 0] - yo[:, 1])
        np.testing.assert_allclose(separation[2:-2], 10.0, rtol=2e-3)

    def test_offset_of_a_circle_has_the_right_radii(self):
        angle = np.linspace(0, 1.5 * np.pi, 900)
        radius = 40.0
        x, y = radius * np.cos(angle), radius * np.sin(angle)
        xo, yo = geometry.offset_polyline(x, y, 4.0)
        r_out = np.hypot(xo[2:-2, 0], yo[2:-2, 0])
        r_in = np.hypot(xo[2:-2, 1], yo[2:-2, 1])
        np.testing.assert_allclose(np.sort([r_out.mean(), r_in.mean()]),
                                   [radius - 4.0, radius + 4.0], rtol=2e-3)

    def test_point_cloud_shape_and_ordering(self):
        line = sine_centerline(n=200, spacing=2.0)
        channel = ChannelConfig(width=20.0, depth=2.0, n_offsets=4)
        topo = bed.compute_bed(line.curvature, line.s, channel,
                               FlipConfig(transverse=False))
        cloud, banks = geometry.build_point_cloud(line.x, line.y, topo, channel,
                                                 collect_banklines=True)
        self.assertEqual(cloud.shape, (200 * 9, 3))
        # The first block is the centerline itself.
        np.testing.assert_allclose(cloud[:200, 0], line.x)
        self.assertIsNotNone(banks)
        polygon = geometry.bankline_polygon(banks)
        np.testing.assert_allclose(polygon[0], polygon[-1])


class TestMesh(unittest.TestCase):
    def test_triangle_count_and_index_range(self):
        n_row, n_off = 20, 3
        triangles = mesh.triangulate(n_row, n_off)
        self.assertEqual(triangles.shape,
                         (2 * (n_row - 1) * 2 * n_off, 3))
        self.assertEqual(triangles.min(), 1)
        self.assertEqual(triangles.max(), n_row * (2 * n_off + 1))

    def test_every_node_is_used(self):
        n_row, n_off = 15, 2
        triangles = mesh.triangulate(n_row, n_off)
        used = np.unique(triangles)
        self.assertEqual(used.size, n_row * (2 * n_off + 1))

    def test_no_degenerate_triangles(self):
        triangles = mesh.triangulate(12, 3)
        for tri in triangles:
            self.assertEqual(len(set(tri.tolist())), 3)

    def test_mesh_files_are_written(self):
        line = sine_centerline(n=60, spacing=5.0)
        channel = ChannelConfig(width=20.0, depth=2.0, n_offsets=3)
        topo = bed.compute_bed(line.curvature, line.s, channel, FlipConfig())
        cloud, _ = geometry.build_point_cloud(line.x, line.y, topo, channel)
        with tempfile.TemporaryDirectory() as tmp:
            written = mesh.write_mesh_files(cloud, line.n_nodes, channel, tmp,
                                            "test")
            self.assertEqual(len(written), 4)
            for path in written:
                self.assertTrue(path.is_file())
                self.assertGreater(path.stat().st_size, 0)
            t3s = [p for p in written if p.suffix == ".t3s"][0]
            self.assertIn(":NodeCount", t3s.read_text().splitlines()[0])


class TestMigration(unittest.TestCase):
    def setUp(self):
        self.channel = ChannelConfig(width=20.0, depth=2.0, n_offsets=5)
        self.rng = np.random.default_rng(0)

    def test_straight_channel_does_not_migrate_without_noise(self):
        x = np.linspace(0, 2000, 400)
        y = np.zeros_like(x)
        s = planform.arc_length(x, y)
        cur, theta = planform.curvature(s, x, y)
        line = Centerline(x=x, y=y, s=s, curvature=cur, theta=theta)
        config = MigrationConfig(enabled=True, ub0=0.0, c0=0.0, e0=1e-7,
                                 dt=86400.0)
        result = migration.migrate(line, cur, cur, self.channel, config,
                                   self.rng)
        np.testing.assert_allclose(result.centerline.y, 0.0, atol=1e-9)
        self.assertAlmostEqual(result.mean_rate, 0.0, places=9)

    def test_curved_channel_migrates(self):
        line = sine_centerline(n=800, spacing=2.0, wavelength=400.0,
                               amplitude=40.0)
        cur_f = planform.filter_curvature(line.curvature)
        cur_l = planform.phase_lag(cur_f, 40)
        config = MigrationConfig(enabled=True, e0=1e-7, dt=86400.0)
        result = migration.migrate(line, cur_f, cur_l, self.channel, config,
                                   self.rng)
        self.assertGreater(result.mean_rate, 0.0)
        self.assertEqual(result.centerline.n_nodes, line.n_nodes)

    def test_lagged_term_dominates_the_local_term(self):
        # This inequality is why meanders grow and translate rather than just
        # oscillating in place, so it is worth pinning down.
        line = sine_centerline(n=400, spacing=2.0)
        config = MigrationConfig(enabled=True, cf0=0.01, fr0=0.1)
        channel = ChannelConfig(width=210.0, depth=5.0, n_offsets=10)
        only_local = migration.near_bank_velocity(
            line, line.curvature, np.zeros_like(line.curvature), channel,
            config, self.rng)
        only_lagged = migration.near_bank_velocity(
            line, np.zeros_like(line.curvature), line.curvature, channel,
            config, self.rng)
        self.assertGreater(np.max(np.abs(only_lagged)),
                           np.max(np.abs(only_local)))

    @staticmethod
    def looping_centerline(theta0=2.6, wavelength=4000.0, ds=20.0, n_bends=4):
        """A sine-generated curve amplitudinous enough to close its own loops.

        The loops touch in the interior of the reach, which is what a real neck
        cutoff looks like.
        """
        s = np.arange(0, n_bends * wavelength, ds)
        angle = theta0 * np.sin(2 * np.pi * s / wavelength)
        x = np.concatenate(([0.0], np.cumsum(ds * np.cos(angle))[:-1]))
        y = np.concatenate(([0.0], np.cumsum(ds * np.sin(angle))[:-1]))
        arc = planform.arc_length(x, y)
        cur, theta = planform.curvature(arc, x, y)
        return Centerline(x=x, y=y, s=arc, curvature=cur, theta=theta)

    def test_neck_cutoff_found_on_a_touching_loop(self):
        # A loop that comes back on itself: two nodes far apart along the
        # channel end up within one channel width of each other.
        line = self.looping_centerline()
        i, j = migration.find_neck_cutoff(line, self.channel,
                                         NeckCutoffConfig())
        self.assertGreaterEqual(i, 0)
        self.assertGreater(j, i)
        self.assertGreater(line.s[j] - line.s[i],
                           2 * self.channel.width)
        self.assertLess(np.hypot(line.x[j] - line.x[i], line.y[j] - line.y[i]),
                        self.channel.width)

    def test_neck_cutoff_end_margin_ignores_boundary_pairs(self):
        # A closed curve: its only touching pair is the two reach ends, which
        # is the boundary artifact the margin exists to reject.
        t = np.linspace(0, 2 * np.pi, 1200)
        x, y = 100 * np.sin(t), 60 * np.sin(2 * t)
        s = planform.arc_length(x, y)
        cur, theta = planform.curvature(s, x, y)
        line = Centerline(x=x, y=y, s=s, curvature=cur, theta=theta)
        self.assertEqual(
            migration.find_neck_cutoff(line, self.channel,
                                       NeckCutoffConfig()),
            (-1, -1))
        # v1.x searched the whole reach and did fire on it.
        i, _ = migration.find_neck_cutoff(
            line, self.channel, NeckCutoffConfig(end_margin_widths=0.0))
        self.assertGreaterEqual(i, 0)

    def test_neck_cutoff_margin_keeps_interior_cutoffs(self):
        line = self.looping_centerline()
        self.assertEqual(
            migration.find_neck_cutoff(line, self.channel,
                                       NeckCutoffConfig()),
            migration.find_neck_cutoff(
                line, self.channel,
                NeckCutoffConfig(end_margin_widths=0.0)))

    def test_end_taper_pins_the_ends_and_frees_the_interior(self):
        line = sine_centerline(n=1200, spacing=2.0)
        config = MigrationConfig(enabled=True, end_taper_widths=2.0)
        weights = migration.end_taper(line, self.channel, config)
        self.assertEqual(weights.size, line.n_nodes)
        self.assertAlmostEqual(weights[0], 0.0)
        self.assertAlmostEqual(weights[-1], 0.0)
        self.assertTrue(np.all((weights >= 0.0) & (weights <= 1.0)))
        self.assertAlmostEqual(np.max(weights), 1.0)
        # Monotone up from the inlet, and symmetric about mid-reach.
        taper = 2.0 * self.channel.width
        inside = line.s <= taper
        self.assertTrue(np.all(np.diff(weights[inside]) >= -1e-12))
        # Symmetric about mid-reach, to within the fixture's slightly
        # non-uniform node spacing.
        np.testing.assert_allclose(weights[1], weights[-2], rtol=1e-2)

    def test_end_taper_off_is_all_ones(self):
        line = sine_centerline(n=200, spacing=2.0)
        weights = migration.end_taper(
            line, self.channel, MigrationConfig(end_taper_widths=0.0))
        np.testing.assert_allclose(weights, 1.0)

    def test_migration_does_not_move_the_end_nodes(self):
        line = sine_centerline(n=600, spacing=5.0)
        config = MigrationConfig(enabled=True, e0=1e-5, dt=86400.0, ub0=3.0)
        result = migration.migrate(line, line.curvature, line.curvature,
                                   self.channel, config, self.rng)
        moved = result.centerline
        self.assertAlmostEqual(moved.x[0], line.x[0], places=9)
        self.assertAlmostEqual(moved.y[0], line.y[0], places=9)
        self.assertAlmostEqual(moved.x[-1], line.x[-1], places=9)
        self.assertAlmostEqual(moved.y[-1], line.y[-1], places=9)
        # The interior did move.
        mid = line.n_nodes // 2
        self.assertNotAlmostEqual(moved.x[mid], line.x[mid], places=9)

    def test_no_neck_cutoff_on_a_gentle_sine(self):
        line = sine_centerline(n=1200, spacing=2.0, wavelength=800.0,
                               amplitude=60.0)
        i, j = migration.find_neck_cutoff(line, self.channel,
                                         NeckCutoffConfig())
        self.assertEqual((i, j), (-1, -1))

    def test_chute_entrances_find_apexes_and_inflections(self):
        s = np.linspace(0, 1000, 2000)
        cur = np.sin(2 * np.pi * s / 200)
        apexes = migration.find_chute_entrances(
            cur, ChuteCutoffConfig(entrance="apex"))
        inflections = migration.find_chute_entrances(
            cur, ChuteCutoffConfig(entrance="inflection"))
        self.assertEqual(apexes.size, inflections.size - 1)
        # Apexes sit between consecutive inflections.
        for k in range(apexes.size):
            self.assertLess(inflections[k], apexes[k])
            self.assertLess(apexes[k], inflections[k + 1])

    def test_chute_valley_angle_bounds(self):
        x = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
        y = np.array([0.0, 1.0, 0.0, 1.0, 0.0])
        # Along the valley axis (which is +x here).
        self.assertAlmostEqual(migration.chute_valley_angle(x, y, 0, 4), 0.0)
        # Perpendicular to it.
        self.assertAlmostEqual(migration.chute_valley_angle(x, y, 0, 1), 45.0)

    def test_chute_cutoff_frequency_is_respected(self):
        line = sine_centerline(n=4000, spacing=1.0, wavelength=400.0,
                               amplitude=60.0)
        cur = planform.filter_curvature(line.curvature)
        channel = ChannelConfig(width=20.0, depth=2.0, n_offsets=5)
        config = ChuteCutoffConfig(enabled=True, frequency=0.5, start_step=0,
                                   max_valley_angle=90.0,
                                   min_length_widths=2.0, min_sinuosity=1.0)
        rng = np.random.default_rng(1)
        hits = sum(1 for _ in range(400)
                   if migration.find_chute_cutoff(line, cur, 10, channel,
                                                  config, rng)[0] >= 0)
        self.assertGreater(hits, 150)
        self.assertLess(hits, 250)

    def test_chute_cutoff_respects_start_step_and_switch(self):
        line = sine_centerline(n=4000, spacing=1.0)
        cur = planform.filter_curvature(line.curvature)
        channel = ChannelConfig(width=20.0, depth=2.0, n_offsets=5)
        rng = np.random.default_rng(2)
        early = ChuteCutoffConfig(enabled=True, frequency=1.0, start_step=100,
                                  min_sinuosity=1.0)
        self.assertEqual(
            migration.find_chute_cutoff(line, cur, 5, channel, early, rng),
            (-1, -1))
        off = ChuteCutoffConfig(enabled=False, frequency=1.0, start_step=0)
        self.assertEqual(
            migration.find_chute_cutoff(line, cur, 500, channel, off, rng),
            (-1, -1))

    def test_chute_cutoff_rejects_a_straight_reach(self):
        # min_sinuosity is what stops a chute being carved across a reach that
        # is already straight, which would gain the flow nothing.
        x = np.linspace(0, 4000, 4000)
        y = 0.01 * np.sin(2 * np.pi * x / 400)      # essentially straight
        s = planform.arc_length(x, y)
        cur, theta = planform.curvature(s, x, y)
        line = Centerline(x=x, y=y, s=s, curvature=cur, theta=theta)
        channel = ChannelConfig(width=20.0, depth=2.0, n_offsets=5)
        config = ChuteCutoffConfig(enabled=True, frequency=1.0, start_step=0,
                                   max_valley_angle=90.0,
                                   min_length_widths=2.0, min_sinuosity=1.2)
        rng = np.random.default_rng(3)
        for _ in range(20):
            self.assertEqual(
                migration.find_chute_cutoff(line, cur, 10, channel, config,
                                            rng),
                (-1, -1))

    def test_carve_cutoff_shortens_and_returns_the_oxbow(self):
        line = sine_centerline(n=500, spacing=1.0)
        cut, event = migration.carve_cutoff(line, 100, 300, "chute", 7)
        self.assertEqual(cut.n_nodes, 500 - (300 - 100 - 1))
        self.assertEqual(event.oxbow_x.size, 199)
        self.assertEqual(event.kind, "chute")
        self.assertEqual(event.step, 7)
        self.assertLess(cut.s[-1], line.s[-1])


class TestModel(unittest.TestCase):
    def _config(self, tmp, **kwargs):
        config = default_config()
        config.channel.ds = 0.1
        config.kinoshita.n_bends = 2
        config.channel.n_offsets = 4
        config.curvature.smoothing_level = 4
        config.output.directory = tmp
        config.output.save_figures = False
        config.output.save_gif = False
        config.output.log_file = ""
        for key, value in kwargs.items():
            section, _, field = key.partition("__")
            if field:
                setattr(getattr(config, section), field, value)
            else:
                setattr(config, section, value)
        return config.validate()

    def test_static_run_writes_expected_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = RiverBedModel(self._config(tmp)).run()
            names = {p.name for p in result.files}
            self.assertIn("kinoshita_topo.xyz", names)
            self.assertIn("kinoshita_boundary.i2s", names)
            self.assertIn("kinoshita_mesh.t3s", names)
            self.assertIn("kinoshita_BC.cli", names)
            self.assertEqual(result.steps_completed, 0)
            self.assertGreater(result.bed.relief, 0)

    def test_migration_run_produces_time_series(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = self._config(tmp, migration__enabled=True,
                                  migration__n_steps=12,
                                  migration__log_every=100,
                                  migration__plot_every=100,
                                  migration__seed=5,
                                  output__save_mesh=False)
            result = RiverBedModel(config).run()
            self.assertEqual(result.sinuosity.size, 12)
            self.assertEqual(result.migration_rate.size, 12)
            self.assertTrue(np.all(result.sinuosity > 1.0))
            self.assertEqual(result.steps_completed, 12)

    def test_seed_makes_a_stochastic_run_reproducible(self):
        with tempfile.TemporaryDirectory() as tmp:
            def once():
                config = self._config(
                    tmp, migration__enabled=True, migration__n_steps=8,
                    migration__ub0=3.0, migration__seed=1234,
                    migration__log_every=100, migration__plot_every=100,
                    output__save_mesh=False, output__save_xyz=False,
                    output__save_bankline=False)
                return RiverBedModel(config).run()
            first, second = once(), once()
            np.testing.assert_allclose(first.centerline.x, second.centerline.x)
            np.testing.assert_allclose(first.sinuosity, second.sinuosity)

    def test_progress_callback_can_stop_a_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = self._config(tmp, migration__enabled=True,
                                  migration__n_steps=50,
                                  migration__log_every=100,
                                  migration__plot_every=100,
                                  output__save_mesh=False)
            result = RiverBedModel(config).run(
                progress=lambda step, total: step < 5)
            self.assertEqual(result.steps_completed, 5)

    def test_two_models_coexist_in_one_process(self):
        # v1.x could not do this: its parameters were module globals baked into
        # Numba kernels at first compile.
        with tempfile.TemporaryDirectory() as tmp:
            narrow = self._config(tmp, channel__width=0.6)
            wide = self._config(tmp, channel__width=6.0)
            a = RiverBedModel(narrow).run()
            b = RiverBedModel(wide).run()
            self.assertAlmostEqual(a.config.channel.width, 0.6)
            self.assertAlmostEqual(b.config.channel.width, 6.0)
            self.assertGreater(b.bed.relief, a.bed.relief)

    def test_centerline_mode_reads_a_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "line.txt"
            s = np.linspace(0, 2000, 400)
            np.savetxt(path, np.column_stack((s, 80 * np.sin(2 * np.pi * s / 500))))
            config = self._config(tmp, mode="centerline",
                                  centerline_file=str(path),
                                  channel__width=40.0, channel__depth=3.0)
            result = RiverBedModel(config).run()
            self.assertEqual(result.config.stem, "line")
            self.assertGreater(result.centerline.sinuosity, 1.0)


class TestCLI(unittest.TestCase):
    def test_init_and_run(self):
        from pyriverbed.cli import main
        with tempfile.TemporaryDirectory() as tmp:
            ini = Path(tmp) / "in.ini"
            self.assertEqual(main(["init", str(ini)]), 0)
            self.assertTrue(ini.is_file())
            self.assertEqual(main(["run", str(ini), "-o", tmp, "--dry-run",
                                   "-q"]), 0)

    def test_init_refuses_to_clobber(self):
        from pyriverbed.cli import main
        with tempfile.TemporaryDirectory() as tmp:
            ini = Path(tmp) / "in.ini"
            main(["init", str(ini)])
            self.assertEqual(main(["init", str(ini)]), 2)
            self.assertEqual(main(["init", str(ini), "--force"]), 0)

    def test_convert_legacy_steering(self):
        from pyriverbed.cli import main
        values = ["jurua.txt", 2, 3, 10.0, 110.0, 0.03125, 0.00520833, 160.0,
                  8.0, 0.0, 0.03, 10, 1, 6.0, 1, 1, 1, 0, 0, 1, 3.0, 0.0, 0.01,
                  0.1, 86400.0, 1e-07, 50, 5000, 100, 24, 20, 1.0, 1, 0.2,
                  1000, 1, 2, 45, 8, 3]
        with tempfile.TemporaryDirectory() as tmp:
            steering = Path(tmp) / "steering.txt"
            steering.write_text("\n".join(str(v) for v in values))
            target = Path(tmp) / "out.ini"
            self.assertEqual(main(["convert", str(steering), str(target)]), 0)
            converted = read_ini(target.read_text())
            self.assertAlmostEqual(converted.channel.width, 160.0)
            self.assertTrue(converted.chute_cutoff.enabled)

    def test_missing_input_file_is_reported(self):
        from pyriverbed.cli import main
        self.assertEqual(main(["run", "/nonexistent/nope.ini", "-q"]), 2)


class TestNotebookAPI(unittest.TestCase):
    def test_configure_shortcuts(self):
        from pyriverbed.notebook import configure
        config = configure(width=1.5, depth=0.3, n_bends=7, chute_cutoff=True,
                           migration=True, chute_frequency=0.25)
        self.assertAlmostEqual(config.channel.width, 1.5)
        self.assertEqual(config.kinoshita.n_bends, 7)
        self.assertTrue(config.chute_cutoff.enabled)
        self.assertAlmostEqual(config.chute_cutoff.frequency, 0.25)

    def test_configure_rejects_unknown_shortcut(self):
        from pyriverbed.notebook import configure
        with self.assertRaises(ConfigError):
            configure(not_a_parameter=1)

    def test_quick_run_and_figures(self):
        import matplotlib
        matplotlib.use("Agg")
        from pyriverbed.notebook import (cross_section, quick_run, show,
                                        show_bed, show_curvature,
                                        show_planform)
        with tempfile.TemporaryDirectory() as tmp:
            result = quick_run(n_bends=2, ds=0.1, n_offsets=4, smoothing=4,
                               output_dir=tmp, save_mesh=False, save_xyz=False,
                               save_figures=False, log_file="", verbose=False)
            for func in (show, show_planform, show_bed, show_curvature,
                         cross_section):
                with self.subTest(func=func.__name__):
                    self.assertIsNotNone(func(result))


class TestArt(unittest.TestCase):
    """The art prints. Every style must render from a real run."""

    @classmethod
    def setUpClass(cls):
        import matplotlib
        matplotlib.use("Agg")
        from pyriverbed.notebook import quick_run
        # e0 * dt must stay a small fraction of a channel width per step:
        # at ~1 width per step the planform blows up, the reach lengthens
        # without bound and the run never finishes.
        cls.result = quick_run(
            n_bends=2, ds=0.1, n_offsets=4, smoothing=4, migration=True,
            n_steps=40, plot_every=10, e0=1e-7, ub0=1.0, seed=3,
            output_dir=tempfile.mkdtemp(), save_mesh=False, save_xyz=False,
            save_figures=False, log_file="", verbose=False)

    def test_every_style_renders(self):
        from pyriverbed import art
        style = art.ArtStyle(dpi=60, paper="a4", bed_grid=90)
        for name in art.STYLES:
            with self.subTest(style=name):
                figure = art.render(self.result, name, style)
                self.assertIsNotNone(figure)
                figure.clf()

    def test_unknown_style_is_rejected(self):
        from pyriverbed import art
        with self.assertRaises(KeyError):
            art.render(self.result, "definitely-not-a-style")

    def test_channel_ribbon_is_a_closed_polygon_of_the_right_width(self):
        from pyriverbed.art import channel_ribbon
        x = np.linspace(0.0, 100.0, 201)
        y = np.zeros_like(x)
        px, py = channel_ribbon(x, y, 3.0)
        self.assertEqual(px.size, py.size)
        self.assertAlmostEqual(px[0], px[-1])
        self.assertAlmostEqual(py[0], py[-1])
        # A straight channel offset by 3 spans exactly 6 across.
        self.assertAlmostEqual(py.max() - py.min(), 6.0, places=6)

    def test_gallery_writes_one_file_per_style(self):
        from pyriverbed import art
        with tempfile.TemporaryDirectory() as tmp:
            written = art.save_gallery(
                self.result, tmp, prefix="t", styles=["minimal", "blueprint"],
                style=art.ArtStyle(dpi=60, paper="a4", bed_grid=90))
            self.assertEqual(len(written), 2)
            for path in written:
                self.assertTrue(path.exists())
                self.assertGreater(path.stat().st_size, 0)

    def test_paper_turns_landscape_for_a_wide_reach(self):
        from pyriverbed.art import ArtStyle
        style = ArtStyle(paper="a4")
        wide, tall = style.figsize(3.0), style.figsize(0.3)
        self.assertGreater(wide[0], wide[1])
        self.assertLess(tall[0], tall[1])


class TestDocumentationMath(unittest.TestCase):
    """The LaTeX in the docs has to survive GitHub's Markdown renderer.

    GitHub renders ``$...$`` and ``$$...$$`` with MathJax, but the Markdown
    parser gets there first, and a few constructs silently come out as literal
    text: inline math straddling a line break, a display block whose ``$$``
    shares a line with the formula, or an unbalanced delimiter. None of those
    raise anything, they just look broken on the rendered page, so they are
    worth a test.
    """

    DOCS = ("THEORY_GUIDE.md", "README.md", "examples/README.md")
    MATH_ENVS = frozenset({"cases", "aligned", "bmatrix", "array", "pmatrix",
                           "matrix", "split", "gather", "align", "vmatrix"})
    #: Macros GitHub's MathJax refuses outright, reporting "The following
    #: macros are not allowed" and rendering the block as raw source.
    DISALLOWED = frozenset({
        "operatorname", "newcommand", "renewcommand", "def", "let", "require",
        "href", "class", "style", "cssId", "label", "tag", "ref", "eqref",
        "includegraphics", "input", "include", "write", "catcode"})

    MATH_FENCE = re.compile(
        r"(?m)^[ \t]*```math[ \t]*\n(.*?)\n[ \t]*```[ \t]*$", re.S)

    @classmethod
    def _split(cls, text):
        """Return the display blocks and the inline spans of *text*.

        Display math is written as a ```math fence rather than as ``$$``:
        GitHub's Markdown parser resolves escapes inside a ``$$`` block, so a
        matrix row separator ``\\`` reaches MathJax as a single backslash and
        the whole block collapses onto one line. A fenced block is passed
        through untouched.
        """
        display = [(text[:m.start()].count("\n") + 1, m.group(1))
                   for m in cls.MATH_FENCE.finditer(text)]
        rest = text
        for _, body in display:
            rest = rest.replace(body, "")
        rest = re.sub(r"```.*?```", "", rest, flags=re.S)
        rest = re.sub(r"`[^`\n]*`", "", rest)
        inline = [(rest[:m.start()].count("\n") + 1, m.group(1))
                  for m in re.finditer(r"(?<!\$)\$([^$\n]+)\$(?!\$)", rest)]
        return display, inline, rest

    def _docs(self):
        root = Path(__file__).resolve().parent.parent
        for name in self.DOCS:
            path = root / name
            if path.exists():
                yield path, path.read_text(encoding="utf-8")

    def test_display_delimiters_are_alone_on_their_lines(self):
        for path, text in self._docs():
            for number, line in enumerate(text.split("\n"), start=1):
                stripped = line.strip()
                if stripped.startswith("$$") and stripped != "$$":
                    self.fail(f"{path.name}:{number}: '$$' shares a line with "
                              f"the formula, so GitHub renders it literally")

    def test_inline_math_never_straddles_a_line_break(self):
        for path, text in self._docs():
            _, _, rest = self._split(text)
            for match in re.finditer(r"(?<!\$)\$([^$]+?)\$(?!\$)", rest):
                if "\n" in match.group(1):
                    line = rest[:match.start()].count("\n") + 1
                    self.fail(f"{path.name}: inline math spans a line break "
                              f"near line {line}")

    def test_inline_math_has_no_padding_spaces(self):
        for path, text in self._docs():
            _, inline, _ = self._split(text)
            for number, body in inline:
                self.assertFalse(
                    body[:1].isspace() or body[-1:].isspace(),
                    f"{path.name}:{number}: '${body}$' is padded with spaces, "
                    f"which stops GitHub treating it as math")

    def test_delimiters_and_environments_balance(self):
        for path, text in self._docs():
            display, inline, _ = self._split(text)
            for number, body in display + inline:
                where = f"{path.name}:{number}"
                self.assertEqual(body.count("{"), body.count("}"),
                                 f"{where}: unbalanced braces")
                self.assertEqual(
                    len(re.findall(r"\\left(?![A-Za-z])", body)),
                    len(re.findall(r"\\right(?![A-Za-z])", body)),
                    f"{where}: \\left without \\right")
                self.assertEqual(body.count("\\lvert"), body.count("\\rvert"),
                                 f"{where}: \\lvert without \\rvert")
                self.assertEqual(body.count("\\lVert"), body.count("\\rVert"),
                                 f"{where}: \\lVert without \\rVert")
                self.assertNotIn("$", body, f"{where}: stray $ inside math")
                for env in re.findall(r"\\begin\{([A-Za-z*]+)\}", body):
                    self.assertIn(env.rstrip("*"), self.MATH_ENVS,
                                  f"{where}: {env!r} is not a math environment")
                    self.assertEqual(body.count(f"\\begin{{{env}}}"),
                                     body.count(f"\\end{{{env}}}"),
                                     f"{where}: {env} begin/end mismatch")

    def test_no_disallowed_macros(self):
        for path, text in self._docs():
            display, inline, _ = self._split(text)
            for number, body in display + inline:
                for macro in set(re.findall(r"\\([A-Za-z]+)", body)):
                    self.assertNotIn(
                        macro, self.DISALLOWED,
                        f"{path.name}:{number}: \\{macro} is on GitHub's "
                        f"disallowed-macro list and will render as raw source")

    def test_display_math_is_fenced_not_dollar_delimited(self):
        # A $$ block has its escapes resolved by the Markdown parser first, so
        # a matrix row separator arrives at MathJax as a single backslash and
        # the block collapses onto one line, often styled as a heading.
        for path, text in self._docs():
            for number, line in enumerate(text.split("\n"), start=1):
                self.assertNotEqual(
                    line.strip(), "$$",
                    f"{path.name}:{number}: use a ```math fence for display "
                    f"math, not $$")

    def test_inline_math_has_no_row_separators(self):
        for path, text in self._docs():
            _, inline, _ = self._split(text)
            for number, body in inline:
                self.assertNotIn(
                    "\\\\", body,
                    f"{path.name}:{number}: '\\\\' inside inline math is eaten "
                    f"by the Markdown parser")

    def test_no_unsupported_delimiter_syntax(self):
        # GitHub accepts $...$ and $$...$$ only, never \( \) or \[ \].
        for path, text in self._docs():
            display, inline, _ = self._split(text)
            for number, body in display + inline:
                self.assertIsNone(
                    re.search(r"(?<!\\)\\[()]", body),
                    f"{path.name}:{number}: \\( \\) is not supported")

    def test_math_contains_nothing_the_markdown_parser_would_eat(self):
        # '<' before a letter reads as an HTML tag and '&word;' as an entity;
        # either would be stripped before MathJax ever saw the formula.
        for path, text in self._docs():
            display, inline, _ = self._split(text)
            for number, body in display + inline:
                where = f"{path.name}:{number}"
                self.assertIsNone(re.search(r"<[A-Za-z/!]", body),
                                  f"{where}: '<' reads as an HTML tag")
                self.assertIsNone(re.search(r"&[A-Za-z]+;", body),
                                  f"{where}: looks like an HTML entity")

    def test_dollar_signs_balance_overall(self):
        for path, text in self._docs():
            _, _, rest = self._split(text)
            self.assertEqual(rest.count("$") % 2, 0,
                             f"{path.name}: odd number of '$' outside display "
                             f"blocks, so one of them opens math by accident")

    def test_theory_guide_documents_the_governing_equations(self):
        # Cheap guard that the guide has not been gutted.
        root = Path(__file__).resolve().parent.parent
        guide = root / "THEORY_GUIDE.md"
        if not guide.exists():
            self.skipTest("THEORY_GUIDE.md not present")
        display, inline, _ = self._split(guide.read_text(encoding="utf-8"))
        self.assertGreaterEqual(len(display), 20)
        self.assertGreaterEqual(len(inline), 100)


if __name__ == "__main__":
    unittest.main(verbosity=2)
