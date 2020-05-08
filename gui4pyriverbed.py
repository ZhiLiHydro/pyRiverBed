# coding: utf-8
"""
Graphical user interface for pyRiverBed
pyRiverBed: Generate Synthetic Riverbed Topography for Meandering Rivers
MIT License
Author: Zhi Li, Univ. of Illinois Urbana-Champaign
Contact: zhil2[at]illinois[dot]edu

"""

import platform
import webbrowser
from shutil import rmtree
from os import remove
from glob import glob
from urllib import request, error
from base64 import encodebytes
import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import pyriverbed


def main():
    """
    Add widgets to the GUI one by one.
    If screen height <1000 pixels (e.g. laptop screens) -> left & right layout
    If screen height >=1000 pixels (e.g. PC screens) -> top & bottom layout
	
    tk and ttk widgets used:
    Canvas, Message, Frame, Label, Scale, Radiobutton, Combobox and Button
    
    Parameters
    ----------
    None
    
    Returns
    -------
    None

    """
    root = tk.Tk()
    root.title('pyRiverBed')
    root.resizable(0, 0)
    
    canvas = tk.Canvas(root)

    if root.winfo_screenheight() < 1051:
        canvas.pack(side='left', padx=20)
    else:
        canvas.pack(side='top', padx=20)

####### MODE SELECTOR 
    msg = tk.Message(
        canvas,
        text='MODE SELECTOR',
        bg='wheat',
        font='System 11 bold',
        width=500,
        relief='raised')
    msg.pack() 

    row = tk.Frame(canvas)
    row.pack(side='top')
    
    tk.Label(
        row,
        text='Mode',
        width=7,
        font='System 10',
        anchor='e').pack(side='left')
    
    modvar = tk.IntVar()
    mod = tk.Radiobutton(
        row,
        text='1=Generate Kinoshita Curve from equation',
        width=43,
        font='System 10',
        anchor='sw',
        variable=modvar,
        value=1)
    mod.select()
    mod.pack() 
    
    mod = tk.Radiobutton(
        row,
        text='2=Read your own centerline from file',
        width=43,
        font='System 10',
        anchor='sw',
        variable=modvar,
        value=2)
    mod.pack() 

    row = tk.Frame(canvas)
    row.pack(side='top')
    
    tk.Label(
        row,
        text='File name',
        width=10,
        font='System 10',
        anchor='c').pack(side='left')
    
    fnamevar = tk.StringVar()
    fname = ttk.Combobox(
        row,
        font='System 10',
        width=18,
        textvariable=fnamevar)
    fname.set('mycenterline.txt')
    fname.pack()
    
    row = tk.Frame(canvas)
    row.pack(side='top')

    if platform.system() == 'Darwin':
        w1, w2 = 15, 120
    else:
        w1, w2 = 20, 135

    tk.Label(
        row,
        text='Level of smoothing ',
        width=w1,
        font='System 10',
        anchor='e').pack(side='left')
    
    smoothvar = tk.IntVar()
    scale = tk.Scale(row, 
                     from_=0, 
                     to=100, 
                     orient='horizontal',
                     width=8,
                     sliderlength=10,
                     font='System 8',
                     tickinterval=25,
                     length=w2,
                     variable=smoothvar)
    scale.set(10)
    scale.pack()

    if platform.system() == 'Darwin':
        row = tk.Frame(canvas)
        row.pack(side='top',pady=5)    

####### KINOSHITA PARAMETERS 
    msg = tk.Message(
        canvas,
        text='KINOSHITA PARAMETERS (ONLY WORK IF MODE=1)',
        bg='light cyan',
        font='System 11 bold',
        width=500,
        relief='raised')
    msg.pack()

    row = tk.Frame(canvas)
    row.pack(side='top')
    
    tk.Label(
        row,
        text='Number of bends (unit: /)',
        width=42,
        font='System 10',
        anchor='w').pack(side='left')
    
    nbendvar = tk.IntVar()
    nbend = ttk.Combobox(
        row,
        font='System 10',
        width=18,
        textvariable=nbendvar)
    nbend['values'] = ('1', '2', '3', '4', '5', '10')
    nbend.set('3')
    nbend.pack()

    row = tk.Frame(canvas)
    row.pack(side='top')
    
    tk.Label(
        row,
        text='Arc wavelength (unit: m)',
        width=42,
        font='System 10',
        anchor='w').pack(side='left')
    
    lamdvar = tk.DoubleVar()
    lamd = ttk.Combobox(
        row,
        font='System 10',
        width=18,
        textvariable=lamdvar)
    lamd['values'] = ('1', '2', '5', '10', '20', '50', '100')
    lamd.set('10')
    lamd.pack()

    row = tk.Frame(canvas)
    row.pack(side='top')
    
    tk.Label(
        row,
        text='Max angular amplitude (unit: deg)',
        width=42,
        font='System 10',
        anchor='w').pack(side='left')
    
    thetavar = tk.DoubleVar()
    thetau = ttk.Combobox(
        row,
        font='System 10',
        width=18,
        textvariable=thetavar)
    thetau['values'] = ('30', '45', '60', '90', '110')
    thetau.set('110')
    thetau.pack()

    row = tk.Frame(canvas)
    row.pack(side='top')
    
    tk.Label(
        row,
        text='Skewness coeff. (unit: /)',
        width=42,
        font='System 10',
        anchor='w').pack(side='left')
    
    skewvar = tk.DoubleVar()
    skew = ttk.Combobox(
        row,
        font='System 10',
        width=18,
        textvariable=skewvar)
    skew['values'] = ('0', '0.005', '0.01', '0.02', '0.03')
    skew.set('0.03125')
    skew.pack()

    row = tk.Frame(canvas)
    row.pack(side='top')
    
    tk.Label(
        row,
        text='Flatness coeff. (unit: /)',
        width=42,
        font='System 10',
        anchor='w').pack(side='left')
    
    flatvar = tk.DoubleVar()
    flat = ttk.Combobox(
        row,
        font='System 10',
        width=18,
        textvariable=flatvar)
    flat['values'] = ('0', '0.005', '0.01', '0.02', '0.03')
    flat.set('0.00520833')
    flat.pack()

####### CHANNEL PARAMETERS
    msg = tk.Message(
        canvas,
        text='CHANNEL PARAMETERS',
        bg='light cyan',
        font='System 11 bold',
        width=500,
        relief='raised')
    msg.pack()

    row = tk.Frame(canvas)
    row.pack(side='top')
    
    tk.Label(
        row,
        text='Channel width (unit: m)',
        width=42,
        font='System 10',
        anchor='w').pack(side='left')
    widvar = tk.DoubleVar()
    wid = ttk.Combobox(
        row,
        font='System 10',
        width=18,
        textvariable=widvar)
    wid.set('0.6')
    wid.pack()

    row = tk.Frame(canvas)
    row.pack(side='top')
    
    tk.Label(
        row,
        text='Water depth (unit: m)',
        width=42,
        font='System 10',
        anchor='w').pack(side='left')
    
    depvar = tk.DoubleVar()
    dep = ttk.Combobox(
        row,
        font='System 10',
        width=18,
        textvariable=depvar)
    dep.set('0.15')
    dep.pack()

    row = tk.Frame(canvas)
    row.pack(side='top')
    
    tk.Label(
        row,
        text='Channel slope (unit: /)',
        width=42,
        font='System 10',
        anchor='w').pack(side='left')
    
    slovar = tk.DoubleVar()
    slo = ttk.Combobox(
        row,
        font='System 10',
        width=18,
        textvariable=slovar)
    slo['values'] = ('0', '1e-3', '5e-4', '1e-4')
    slo.set('0')
    slo.pack()
    
    
    row = tk.Frame(canvas)
    row.pack(side='top')
    
    tk.Label(
        row,
        text='Transverse slope corrector (unit: /)',
        width=42,
        font='System 10',
        anchor='w').pack(side='left')
    
    stvar = tk.DoubleVar()
    st = ttk.Combobox(
        row,
        font='System 10',
        width=18,
        textvariable=stvar)
    st['values'] = ('0.5', '0.6', '0.7', '0.8', '0.9', '1')
    st.set('1')
    st.pack()

    row = tk.Frame(canvas)
    row.pack(side='top')
    
    tk.Label(
        row,
        text='Streamwise resolution (unit: m)',
        width=42,
        font='System 10',
        anchor='w').pack(side='left')
    
    dsvar = tk.DoubleVar()
    ds = ttk.Combobox(
        row,
        font='System 10',
        width=18,
        textvariable=dsvar)
    ds.set('0.03')
    ds.pack()
    
    row = tk.Frame(canvas)
    row.pack(side='top')
    
    tk.Label(
        row,
        text='(Grid size of centerline, ONLY WORK IF MODE=1)',
        width=55,
        font='System 10',
        anchor='w').pack(side='left')

    row = tk.Frame(canvas)
    row.pack(side='top')
    
    tk.Label(
        row,
        text='Transverse resolution (unit: /)',
        width=42,
        font='System 10',
        anchor='w').pack(side='left')
    
    noffsvar = tk.IntVar()
    noffs = ttk.Combobox(
        row,
        font='System 10',
        width=18,
        textvariable=noffsvar)
    noffs['values'] = ('5', '8', '10', '20')
    noffs.set('10')
    noffs.pack()
    
    row = tk.Frame(canvas)
    row.pack(side='top')
    
    tk.Label(
        row,
        text='(Number of polyline offsets at each side of centerline)',
        width=55,
        font='System 10',
        anchor='w').pack(side='left')

    canvas = tk.Canvas(root)
    if root.winfo_screenheight() < 1000:
        canvas.pack(side='left', padx=20)
    else:
        canvas.pack(side='top', padx=20)

####### PHASE LAG 
    msg = tk.Message(
        canvas,
        text='CURVATURE PHASE LAG PARAMETERS',
        bg='light cyan',
        font='System 11 bold',
        width=500,
        relief='raised')
    msg.pack()

    row = tk.Frame(canvas)
    row.pack(side='top')
    
    tk.Label(
        row,
        text='Curvature phase lag',
        width=19,
        font='System 10',
        anchor='c').pack(side='left')
    
    if platform.system() == 'Darwin':
        w = 8
    else:
        w = 4
        
    lagmodvar = tk.IntVar()
    lagmod = tk.Radiobutton(
        row,
        text='Off',
        width=w,
        font='System 10',
        anchor='sw',
        variable=lagmodvar,
        value=0)
    lagmod.pack(side='left')
    
    lagmod = tk.Radiobutton(
        row,
        text='On',
        width=w,
        font='System 10',
        anchor='sw',
        variable=lagmodvar,
        value=1)
    lagmod.select()
    lagmod.pack(side='left')   
    
    tk.Label(
        row,
        text='Lag strength (unit: /)',
        width=19,
        font='System 10',
        anchor='w').pack(side='left')
    
    lagstrthvar = tk.DoubleVar()
    lagstrth = ttk.Combobox(
        row,
        font='System 10',
        width=3,
        textvariable=lagstrthvar)
    lagstrth['values'] = ('2', '3', '4', '5', '6')
    lagstrth.set('4')
    lagstrth.pack(side='left')

####### OUTPUT
    msg = tk.Message(
        canvas,
        text='FILE OUTPUT',
        bg='light cyan',
        font='System 11 bold',
        width=500,
        relief='raised')
    msg.pack()

    row = tk.Frame(canvas)
    row.pack(side='top')
    
    if platform.system() == 'Darwin':
        w = 10
    else:
        w = 7
    
    tk.Label(
        row,
        text='Save riverbed topography file (.xyz)?',
        width=42,
        font='System 10',
        anchor='w').pack(side='left')

    xyzvar = tk.IntVar()
    xyz = tk.Radiobutton(
        row,
        text='Yes',
        width=w,
        font='System 10',
        anchor='sw',
        variable=xyzvar,
        value=1)
    xyz.pack(side='left') 
    xyz.select()
    xyz = tk.Radiobutton(
        row,
        text='No',
        width=w,
        font='System 10',
        anchor='sw',
        variable=xyzvar,
        value=0)
    xyz.pack(side='left') 

    row = tk.Frame(canvas)
    row.pack(side='top')
    
    tk.Label(
        row,
        text='Save river bankline file (.i2s)?',
        width=42,
        font='System 10',
        anchor='w').pack(side='left')
    
    bdvar = tk.IntVar()
    bd = tk.Radiobutton(
        row,
        text='Yes',
        width=w,
        font='System 10',
        anchor='sw',
        variable=bdvar,
        value=1)
    bd.pack(side='left') 
    bd.select()
    bd = tk.Radiobutton(
        row,
        text='No',
        width=w,
        font='System 10',
        anchor='sw',
        variable=bdvar,
        value=0)
    bd.pack(side='left') 
    
    row = tk.Frame(canvas)
    row.pack(side='top')
    
    tk.Label(
        row,
        text='Save FEM mesh & BC files (.t3s, .dat, .cli, .bc2)?',
        width=42,
        font='System 10',
        anchor='w').pack(side='left')
    
    meshvar = tk.IntVar()
    mesh = tk.Radiobutton(
        row,
        text='Yes',
        width=w,
        font='System 10',
        anchor='sw',
        variable=meshvar,
        value=1)
    mesh.pack(side='left') 
    mesh.select()
    mesh = tk.Radiobutton(
        row,
        text='No',
        width=w,
        font='System 10',
        anchor='sw',
        variable=meshvar,
        value=0)
    mesh.pack(side='left') 
    
####### FLIP 
    msg = tk.Message(
        canvas,
        text='FLIP',
        bg='light cyan',
        font='System 11 bold',
        width=500,
        relief='raised')
    msg.pack()

    row = tk.Frame(canvas)
    row.pack(side='top')
    
    tk.Label(
        row,
        text='Flip in streamwise direction?',
        width=42,
        font='System 10',
        anchor='w').pack(side='left')

    flpstrmvar = tk.IntVar()
    flpstrm = tk.Radiobutton(
        row,
        text='Yes',
        width=w,
        font='System 10',
        anchor='sw',
        variable=flpstrmvar,
        value=1)
    flpstrm.pack(side='left') 
    flpstrm = tk.Radiobutton(
        row,
        text='No',
        width=w,
        font='System 10',
        anchor='sw',
        variable=flpstrmvar,
        value=0)
    flpstrm.pack(side='left')
    flpstrm.select()
    
    row = tk.Frame(canvas)
    row.pack(side='top')
    
    tk.Label(
        row,
        text='Flip in transverse direction?',
        width=42,
        font='System 10',
        anchor='w').pack(side='left')

    flptranvar = tk.IntVar()
    flptran = tk.Radiobutton(
        row,
        text='Yes',
        width=w,
        font='System 10',
        anchor='sw',
        variable=flptranvar,
        value=1)
    flptran.pack(side='left') 
    flptran.select()
    flptran = tk.Radiobutton(
        row,
        text='No',
        width=w,
        font='System 10',
        anchor='sw',
        variable=flptranvar,
        value=0)
    flptran.pack(side='left')

####### MIGRATION 
    msg = tk.Message(
        canvas,
        text='MEANDER MIGRATION PARAMETERS',
        bg='light cyan',
        font='System 11 bold',
        width=500,
        relief='raised')
    msg.pack()

    row = tk.Frame(canvas)
    row.pack(side='top')
    
    if platform.system() == 'Darwin':
        w = 8
    else:
        w = 4
        
    tk.Label(
        row,
        text='Meander migration',
        width=17,
        font='System 10',
        anchor='c').pack(side='left')
    migvar = tk.IntVar()
    mig = tk.Radiobutton(
        row,
        text='Off',
        width=w,
        font='System 10',
        anchor='sw',
        variable=migvar,
        value=0)
    mig.select()
    mig.pack(side='left')
    
    mig = tk.Radiobutton(
        row,
        text='On',
        width=w,
        font='System 10',
        anchor='sw',
        variable=migvar,
        value=1)
    mig.pack(side='left')

    row = tk.Frame(canvas)
    row.pack(side='top')
    
    tk.Label(
        row,
        text='Ub0 (unit: /)',
        width=15,
        font='System 10',
        anchor='w').pack(side='left')
    ub0var = tk.DoubleVar()
    ub0 = ttk.Combobox(
        row,
        font='System 10',
        width=12,
        textvariable=ub0var)
    ub0.set('0')
    ub0.pack(side='left')

    c0var = tk.DoubleVar()
    c0 = ttk.Combobox(
        row,
        font='System 10',
        width=12,
        textvariable=c0var)
    c0.set('0')
    c0.pack(side='right')
    tk.Label(
        row,
        text='        C0 (unit: /)',
        width=15,
        font='System 10',
        anchor='w').pack(side='right')

    row = tk.Frame(canvas)
    row.pack(side='top')
    
    tk.Label(
        row,
        text='Cf0 (unit: /)',
        width=15,
        font='System 10',
        anchor='w').pack(side='left')
    cf0var = tk.DoubleVar()
    cf0 = ttk.Combobox(
        row,
        font='System 10',
        width=12,
        textvariable=cf0var)
    cf0['values'] = ('0.01', '0.02', '0.03')
    cf0.set('0.01')
    cf0.pack(side='left')

    fr0var = tk.DoubleVar()
    fr0 = ttk.Combobox(
        row,
        font='System 10',
        width=12,
        textvariable=fr0var)
    fr0['values'] = ('0.1', '0.2', '0.3', '0.4', '0.5')
    fr0.set('0.1')
    fr0.pack(side='right')
    tk.Label(
        row,
        text='        Fr0 (unit: /)',
        width=15,
        font='System 10',
        anchor='w').pack(side='right')

    row = tk.Frame(canvas)
    row.pack(side='top')
    
    tk.Label(
        row,
        text='\u0394t (unit: s)',
        width=15,
        font='System 10',
        anchor='w').pack(side='left')
    dtvar = tk.DoubleVar()
    dt = ttk.Combobox(
        row,
        font='System 10',
        width=12,
        textvariable=dtvar)
    dt['values'] = ('60', '3600', '86400')
    dt.set('86400')
    dt.pack(side='left')

    e0var = tk.DoubleVar()
    e0 = ttk.Combobox(
        row,
        font='System 10',
        width=12,
        textvariable=e0var)
    e0['values'] = ('1e-5', '1e-6', '1e-7', '1e-8')
    e0.set('1e-7')
    e0.pack(side='right')
    tk.Label(
        row,
        text='        E0 (unit: /)',
        width=15,
        font='System 10',
        anchor='w').pack(side='right')

    row = tk.Frame(canvas)
    row.pack(side='top')
    
    tk.Label(
        row,
        text='Listing printout period',
        width=20,
        font='System 10',
        anchor='w').pack(side='left')
    Lprintvar = tk.IntVar()
    Lprint = ttk.Combobox(
        row,
        font='System 10',
        width=7,
        textvariable=Lprintvar)
    Lprint.set('50')
    Lprint.pack(side='left')


    tstepsvar = tk.IntVar()
    tsteps = ttk.Combobox(
        row,
        font='System 10',
        width=7,
        textvariable=tstepsvar)
    tsteps.set('10000')
    tsteps.pack(side='right')
    
    tk.Label(
        row,
        text='        # of time steps',
        width=20,
        font='System 10',
        anchor='w').pack(side='right')

    row = tk.Frame(canvas)
    row.pack(side='top')

    tk.Label(
        row,
        text='Graphic printout period',
        width=20,
        font='System 10',
        anchor='w').pack(side='left')
    gprintvar = tk.IntVar()
    gprint = ttk.Combobox(
        row,
        font='System 10',
        width=7,
        textvariable=gprintvar)
    gprint.set('100')
    gprint.pack(side='left')

    fpsvar = tk.IntVar()
    fps = ttk.Combobox(
        row,
        font='System 10',
        width=7,
        textvariable=fpsvar)
    fps['values'] = ('20', '24', '30', '60')
    fps.set('24')
    fps.pack(side='right')
    tk.Label(
        row,
        text='        FPS of GIF',
        width=20,
        font='System 10',
        anchor='w').pack(side='right')

######## BUTTONS 
    row = tk.Frame(canvas)
    row.pack(side='top',pady=3)
    
    row = tk.Frame(canvas)
    row.pack(side='top')

    def gen():
        entries = [modvar.get(), nbendvar.get(), lamdvar.get(),
                   thetavar.get(), skewvar.get(), flatvar.get(),
                   widvar.get(), depvar.get(), slovar.get(),
                   dsvar.get(), noffsvar.get(), lagmodvar.get(),
                   lagstrthvar.get(), xyzvar.get(), 
                   bdvar.get(), meshvar.get(), flpstrmvar.get(), 
                   flptranvar.get(), migvar.get(), ub0var.get(),
                   c0var.get(), cf0var.get(), fr0var.get(),
                   dtvar.get(), e0var.get(), Lprintvar.get(),
                   tstepsvar.get(), gprintvar.get(), fpsvar.get(),
                   smoothvar.get(), stvar.get()]
        with open('steering.txt', 'w') as f:
            entries_str = str(entries)
            f.write(fnamevar.get() + '\n')
            f.write(entries_str.replace(', ', '\n')[1:-1] + '\n')
        messagebox.showinfo(
            message='The steering file \'steering.txt\' has been generated',
            icon='info')

    tk.Button(
            row, 
            text='Generate steering file', 
            font='System 11 bold', 
            command=gen,
            padx=6,
            pady=1).pack(side='left',padx=6,pady=1)

    def run(): 
        if messagebox.showinfo(
            message='pyRiverBed will run',
            icon='info') == 'ok':
            pyriverbed.main()
            
    tk.Button(
        row,
        text='Run pyRiverBed',
        font='System 11 bold',
        command=run,
        padx=6,
        pady=1).pack(side='left',padx=6,pady=1)


    def clean(): 
        if messagebox.askyesno(
            message='Existing *.jpg, *.pdf, *,gif, *.xyz, *.i2s, *.t3s, *.dat, *.cli and *.bc2 files will be deleted',
            icon='warning'):
            extList = ['*.jpg', '*.pdf', '*.gif', '*.xyz', '*.i2s', '*.t3s', '*.dat', '*.cli', '*.bc2']
            n_file_removed = 0
            file_names = ''
            for ext in extList:
                for f in glob(ext):
                    remove(f)
                    n_file_removed += 1
                    file_names += f + '\n'
            if n_file_removed != 0:
                messagebox.showinfo(
                    message = str(n_file_removed) + ' files deleted:\n' + file_names,
                    icon='info')
            else:
                messagebox.showinfo(
                    message='Nothing to delete',
                    icon='info')
            
    tk.Button(
        row,
        text='Clean',
        font='System 11 bold',
        command=clean,
        padx=10,
        pady=1).pack(side='left',padx=6,pady=1)

    def q(): 
        if messagebox.askyesno(
                message='Quit?',
                icon='warning'):
            root.destroy()
            
    tk.Button(
        row, 
        text='Quit', 
        font='System 11 bold', 
        command=q,
        padx=10,
        pady=1).pack(side='left',padx=6,pady=1)

    row = tk.Frame(canvas)
    row.pack(side='top')
 
    def open_github(): 
        webbrowser.open('https://github.com/ZhiLiHydro')

    tk.Label(
        row,
        text='visit',
        width=5,
        font='System 10',
        anchor='e').pack(side='left')

    urlok = True
    try:
        u = request.urlopen('https://github.githubassets.com/images/modules/logos_page/GitHub-Logo.png')
    except error.URLError:
        urlok = False
    if urlok:
        icon = u.read()
        u.close()
        icon = tk.PhotoImage(data=encodebytes(icon)).subsample(15, 15)
        tk.Button(
            row,
            font='System 11 bold', 
            command=open_github,
            width=60,
            height=25,
            image=icon).pack(side='left')
    else:      
        tk.Button(
            row,
            text='GitHub',
            font='System 11 bold', 
            command=open_github,
            padx=15,
            pady=1).pack(side='left',padx=3,pady=1)

    tk.Label(
        row,
        text='for README',
        width=10,
        font='System 10',
        anchor='w').pack(side='left')
    
    root.mainloop()
    
    try:
        rmtree('__pycache__')
    except OSError:
        pass

if __name__ == '__main__':
    main()
    

