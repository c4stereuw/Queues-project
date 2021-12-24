from tkinter import *
from PIL import ImageTk,Image
from tkinter import messagebox
from tkinter import filedialog
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

root=Tk()
root.title('Queues Moduels')
root.geometry("950x521")
filename=PhotoImage(file=Path(__file__).with_name('Qe.png'))
background_label =Label(root,image=filename)
background_label.place(x=0,y=0)


clicked1=StringVar()
clicked1.set("Select type")

clicked2=StringVar()
clicked2.set("Select type")


def deter() :
    global drop1
    drop1=OptionMenu(root,clicked1,"D/D/1/K-1","D/D/1/K-1(M)")
    drop1.grid(row=2 , column=1,padx=15,pady=5)
    drop1.config(bg="#081121",fg="white")
    check = False
    try: drop2
    except NameError: check = True
    if  not check: 
        OptionMenu.destroy(drop2)


def stoch() :
    global drop2
    drop2=OptionMenu(root,clicked2,"M/M/1","M/M/1/K","M/M/C","M/M/C/K")
    drop2.grid(row=2 , column=2,padx=15,pady=5)
    drop2.config(bg="#081121",fg="white")
    check = False
    try: drop1
    except NameError: check = True
    if  not check: 
        OptionMenu.destroy(drop1)

Type_of_Queue = IntVar()

op1= Radiobutton(root,text="Deterministic queue",variable=Type_of_Queue,value=1,padx=1,pady=1,command=deter,bg="#1c1c1b",fg="#ffffff",borderwidth=0,font=38)
op2= Radiobutton(root,text="Stochastic queue",variable=Type_of_Queue,value=2,padx=1,pady=1,command=stoch,bg="#1c1c1b",fg="#ffffff",borderwidth=0,font=38)
op1.grid(row=1 , column=1,padx=165,pady=50)
op2.grid(row=1 , column=2,padx=165,pady=50)

frame0 = LabelFrame(root, text="D/D/1/K-1", padx = 50, pady = 30,bg="#202024",fg="white")
frame1 = LabelFrame(root, text="D/D/1/K-1(M)",padx=50,pady=30,bg="#202024",fg="white")
frame2 = LabelFrame(root, text="M/M/1",padx=50,pady=30,bg="#202024",fg="white")
frame3 = LabelFrame(root, text="M/M/C",padx=50,pady=30,bg="#202024",fg="white")
frame4 = LabelFrame(root, text="M/M/1/K",padx=50,pady=30,bg="#202024",fg="white")
frame5 = LabelFrame(root, text="M/M/C/K",padx=50,pady=30,bg="#202024",fg="white")


def frames() :
    if(Type_of_Queue.get()==1) :
        frame2.grid_forget()
        frame3.grid_forget()
        frame4.grid_forget()
        frame5.grid_forget()
        if(clicked1.get()=="D/D/1/K-1(M)") :
            ti=IntVar()
            nt=IntVar()
            wq=IntVar()
            def DD1KM():
                global lamda
                lamda=1/float(frame1e1.get())
                global mu
                mu=1/float(frame1e2.get())
                k=1+int(frame1e3.get())
                global M
                M=int(frame1e4.get())
                if(mu-lamda>0):
                    ti0=int(M/(mu-lamda))
                    while(M==int((ti0*mu)-int(ti0*lamda))):
                        ti.set(ti0)
                        ti0-=int(frame1e2.get())
                else:
                    ti.set(0)
                T=int(frame1e5.get())
                if(T<ti.get()) :
                    nt0=M+int(T*lamda)-int(T*mu)
                    nt.set(nt0)
                elif (ti.get()==0):
                    nt.set(M)
                else :
                    if(T%1/lamda==0):
                        nt.set(1)
                    else :
                        nt.set(0)
                N=int(frame1e6.get())
                if(N==0) :
                    wq0=(M-1)/(2*mu)
                    wq.set(wq0)
                elif(N<=int(lamda*ti.get())):
                    wq0=(M-1+N)/mu-N/lamda
                    wq.set(wq0)
                else :
                    wq.set(0)
            def graph() :
                arrx =[]
                arry =[]
                for i in range (0,25*int(frame1e2.get()),1):
                    x = i
                    arrx.append(x)
                    if(ti.get()==0):
                        arry.append(M)
                    elif(x<ti.get()) :
                        y=M+int(x*lamda)-int(x*mu)
                        arry.append(y)
                    else :
                        if(x%(1/lamda)==0):
                            y=1
                            arry.append(y)
                        else :
                            y=0
                            arry.append(y)
                plt.step(arrx,arry,data=None,where='post')
                plt.show()

            frame0.grid_forget()
            frame1.grid(row = 4, column=1,columnspan=4,padx=160)
            frame1e1=Entry(frame1,width=30,borderwidth=4,bg="#202024",fg="white")
            frame1e1.grid(row=1,column=2,columnspan=1)
            frame1e2=Entry(frame1,width=30,borderwidth=4,bg="#202024",fg="white")
            frame1e2.grid(row=2,column=2,columnspan=1)
            frame1e3=Entry(frame1,width=30,borderwidth=4,bg="#202024",fg="white")
            frame1e3.grid(row=3,column=2,columnspan=1)
            frame1e4=Entry(frame1,width=30,borderwidth=4,bg="#202024",fg="white")
            frame1e4.grid(row=1,column=4,columnspan=1)
            frame1e5=Entry(frame1,width=30,borderwidth=4,bg="#202024",fg="white")
            frame1e5.grid(row=2,column=4,columnspan=1)
            frame1e6=Entry(frame1,width=30,borderwidth=4,bg="#202024",fg="white")
            frame1e6.grid(row=3,column=4,columnspan=1)
            e1=Label(frame1,text="1/λ",bg="#202024",fg="white").grid(row=1,column=1)
            e2=Label(frame1,text="1/μ",bg="#202024",fg="white").grid(row=2,column=1)
            e3=Label(frame1,text="K-1",bg="#202024",fg="white").grid(row=3,column=1)
            e4=Label(frame1,text="M",bg="#202024",fg="white").grid(row=1,column=3,sticky=E)
            e5=Label(frame1,text="T",bg="#202024",fg="white").grid(row=2,column=3,sticky=E)
            e6=Label(frame1,text="N",bg="#202024",fg="white").grid(row=3,column=3,sticky=E)
            btn0=Button(frame1,text="Calculate",padx=30,command=DD1KM)
            btn0.grid(row = 4, column=2,columnspan=1,pady=6)
            btn1=Button(frame1,text="Graph",padx=40,command=graph)
            btn1.grid(row = 4, column=4,columnspan=1,pady=6)
            e7=Label(frame1,text="ti",bg="#202024",fg="white").grid(row=5,column=1,sticky=E)
            e8=Label(frame1,text="n(t)",bg="#202024",fg="white").grid(row=5,column=3,sticky=E)
            e9=Label(frame1,text="Wq(n)",bg="#202024",fg="white").grid(row=6,column=2,sticky=E)
            e10=Label(frame1,textvariable=ti,padx=80,bg="#327ba8").grid(row=5,column=2,pady=5)
            e11=Label(frame1,textvariable=nt,padx=80,bg="#327ba8").grid(row=5,column=4,pady=5)
            e12=Label(frame1,textvariable=wq,padx=80,bg="#327ba8").grid(row=6,column=3,pady=5)
            
            
            
        elif(clicked1.get()=="D/D/1/K-1") :
            ti=IntVar()
            nt=IntVar()
            wq=IntVar()
            def DD1K():
                global lamda
                global mu
                lamda=1/float(frame0e1.get())
                mu=1/float(frame0e2.get())
                k=1+int(frame0e3.get())
                global capacity1 
                capacity1 =k
                ti0=int((k-mu/lamda)/(lamda-mu)+.01)
                while(k==int((ti0*lamda)-int((mu*ti0-mu/lamda)+.000001))):
                    ti.set(ti0)
                    ti0-=int(frame0e1.get())
                T=int(frame0e5.get())
                if(T<int(frame0e1.get())) :
                    nt.set(0)
                elif(T>=int(frame0e1.get()) and T<ti.get()):
                    nt0=int(T*lamda)
                    nt1=int(mu*(T-1/lamda))
                    nt.set(nt0-nt1)
                else :
                    if(int(frame0e2.get())%int(frame0e1.get())==0):
                        nt.set(k-1)
                    else:
                        T1=T-ti.get()-(1/mu-1/lamda)
                        T2=T1-(1/mu-1/lamda)
                        indeed=int(frame0e2.get()) * (int(frame0e2.get())-int(frame0e1.get()))
                        check = False
                        i=0
                        while(check==False) :
                            if(T1-indeed*i>=0 and T2-indeed*i<0) :
                                check = True
                            elif(T1-indeed*i>=0 and T2-indeed*i>=0) :
                                i+=1
                            else :
                                break
                        if(check) :
                            nt.set(k-2)
                        else :
                            nt.set(k-1)

                N=int(frame0e6.get())
                if(N==0) :
                    wq.set(0)
                elif(N<lamda*ti.get()):
                    wq0=int(frame0e2.get())-int(frame0e1.get())
                    wq.set(wq0*(N-1))
                else :
                    n=lamda*ti.get()
                    wq1=int(frame0e2.get())-int(frame0e1.get())
                    if(int(frame0e2.get())%int(frame0e1.get())==0):
                        wq.set(wq1*(n-2))
                    else:
                        wq.set(str(wq1*(n-2))+ " or "+str(wq1*(n-3)))
            def graph():
                arrx = []
                arry = []
                for i in range (0,25*int(frame0e1.get()),1) :
                    x = i
                    arrx.append(x)
                    if(x<int(frame0e1.get())) :
                        y=0
                        arry.append(y)
                    elif(x>=int(frame0e1.get()) and x<ti.get()):
                        nt0=int(x*lamda)
                        nt1=int(mu*(x-1/lamda))
                        y=nt0-nt1
                        arry.append(y)
                    else :
                        if(int(frame0e2.get())%int(frame0e1.get())==0):
                            y=capacity1-1
                            arry.append(y)
                        else:
                            T1=x-ti.get()-(1/mu-1/lamda)
                            T2=T1-(1/mu-1/lamda)
                            indeed=int(frame0e2.get()) * (int(frame0e2.get())-int(frame0e1.get()))
                            check = False
                            j=0
                            while(check==False) :
                                if(T1-indeed*j>=0 and T2-indeed*j<0) :
                                    check = True
                                elif(T1-indeed*j>=0 and T2-indeed*j>=0) :
                                    j+=1
                                else :
                                    break
                            if(check) :
                                y=capacity1-2
                                arry.append(y)
                            else :
                                y=capacity1-1
                                arry.append(y)
                plt.step(arrx,arry,data=None,where='post')
                plt.show()

            frame1.grid_forget()
            frame0.grid(row = 4, column=1,columnspan=4,padx=160)
            frame0e1=Entry(frame0,width=30,borderwidth=4,bg="#202024",fg="white")
            frame0e1.grid(row=1,column=2,columnspan=1)
            frame0e2=Entry(frame0,width=30,borderwidth=4,bg="#202024",fg="white")
            frame0e2.grid(row=2,column=2,columnspan=1)
            frame0e3=Entry(frame0,width=30,borderwidth=4,bg="#202024",fg="white")
            frame0e3.grid(row=3,column=2,columnspan=1)
            frame0e4=Entry(frame0,width=30,borderwidth=4,state=DISABLED,disabledbackground="#403f4a",fg="white")
            frame0e4.grid(row=1,column=4,columnspan=1)
            frame0e5=Entry(frame0,width=30,borderwidth=4,bg="#202024",fg="white")
            frame0e5.grid(row=2,column=4,columnspan=1)
            frame0e6=Entry(frame0,width=30,borderwidth=4,bg="#202024",fg="white")
            frame0e6.grid(row=3,column=4,columnspan=1)
            e1=Label(frame0,text="1/λ",bg="#202024",fg="white").grid(row=1,column=1)
            e2=Label(frame0,text="1/μ",bg="#202024",fg="white").grid(row=2,column=1)
            e3=Label(frame0,text="K-1",bg="#202024",fg="white").grid(row=3,column=1)
            e4=Label(frame0,text="M",bg="#202024",fg="white").grid(row=1,column=3,sticky=E)
            e5=Label(frame0,text="T",bg="#202024",fg="white").grid(row=2,column=3,sticky=E)
            e6=Label(frame0,text="N",bg="#202024",fg="white").grid(row=3,column=3,sticky=E)
            btn0=Button(frame0,text="Calculate",padx=30,command=DD1K)
            btn0.grid(row = 4, column=2,columnspan=1,pady=6)
            btn1=Button(frame0,text="Graph",padx=40,command=graph)
            btn1.grid(row = 4, column=4,columnspan=1,pady=6)
            e7=Label(frame0,text="ti",bg="#202024",fg="white").grid(row=5,column=1,sticky=E)
            e8=Label(frame0,text="n(t)",bg="#202024",fg="white").grid(row=5,column=3,sticky=E)
            e9=Label(frame0,text="Wq(n)",bg="#202024",fg="white").grid(row=6,column=2,sticky=E)
            e10=Label(frame0,textvariable=ti,padx=80,bg="#327ba8").grid(row=5,column=2,pady=5)
            e11=Label(frame0,textvariable=nt,padx=80,bg="#327ba8").grid(row=5,column=4,pady=5)
            e12=Label(frame0,textvariable=wq,padx=80,bg="#327ba8").grid(row=6,column=3,pady=5)
           
                


    else :
        frame0.grid_forget()
        frame1.grid_forget()

        if(clicked2.get()=="M/M/1") :
            L = IntVar()
            Lq= IntVar()
            W = IntVar()
            Wq= IntVar()
            def MM1() :
                global lamda
                global mu
                lamda=float(frame2e1.get())
                mu=float(frame2e2.get())
                w0 = 1/(float(mu)-float(lamda))
                W.set(round(w0,3))
                Wq.set(round((w0*lamda/mu),3))
                L.set(round((w0*lamda),3))
                Lq.set(round((L.get()*lamda/mu),3))

            frame3.grid_forget()
            frame4.grid_forget()
            frame5.grid_forget()
            frame2.grid(row = 4, column=1,columnspan=4,padx=160)
            frame2e1=Entry(frame2,width=30,borderwidth=4,bg="#202024",fg="white")
            frame2e1.grid(row=1,column=2,columnspan=1)
            frame2e2=Entry(frame2,width=30,borderwidth=4,bg="#202024",fg="white")
            frame2e2.grid(row=1,column=4,columnspan=1)
            frame2e3=Entry(frame2,width=30,borderwidth=4,state=DISABLED,disabledbackground="#403f4a",fg="white")
            frame2e3.grid(row=2,column=3,columnspan=1,pady=10)
            e1=Label(frame2,text="λ",bg="#202024",fg="white").grid(row=1,column=1)
            e2=Label(frame2,text="μ",bg="#202024",fg="white").grid(row=1,column=3,sticky=E)
            e3=Label(frame2,text="c",bg="#202024",fg="white").grid(row=2,column=2,sticky=E)
            btn0=Button(frame2,text="Calculate",padx=30,command=MM1)
            btn0.grid(row = 3, column=3,columnspan=1,pady=6)
            e4=Label(frame2,text="L",bg="#202024",fg="white").grid(row=4,column=1,sticky=E)
            e5=Label(frame2,text="Lq",bg="#202024",fg="white").grid(row=4,column=3,sticky=E)
            e6=Label(frame2,text="W",bg="#202024",fg="white").grid(row=5,column=1,sticky=E)
            e7=Label(frame2,text="Wq",bg="#202024",fg="white").grid(row=5,column=3,sticky=E)
            e8=Label(frame2,textvariable=L,padx=80,bg="#327ba8").grid(row=4,column=2,pady=5)
            e9=Label(frame2,textvariable=Lq,padx=80,bg="#327ba8").grid(row=4,column=4,pady=5)
            e10=Label(frame2,textvariable=W,padx=80,bg="#327ba8").grid(row=5,column=2,pady=5)
            e11=Label(frame2,textvariable=Wq,padx=80,bg="#327ba8").grid(row=5,column=4,pady=5)

        elif(clicked2.get()=="M/M/C") :
            L = IntVar()
            Lq= IntVar()
            W = IntVar()
            Wq= IntVar()
            def MMC() :
                global lamda
                global mu
                global C
                lamda=float(frame3e1.get())
                mu=float(frame3e2.get())
                r = lamda/mu
                C=int(frame3e3.get())
                roh=r/C
                factC=1
                for fac in range(1,C+1):
                    factC=factC*fac
                p00=0
                if(roh<1):
                    for n in range(0,C):
                        factn=1
                        for f in range(1,n+1):
                            factn=factn*f
                        p00+=(r**n/factn)
                    p00+=(C*(r**C)/factC*(C-r))
                else :
                    for n in range(0,C):
                        factn=1
                        for f in range(1,n+1):
                            factn=factn*f
                        p00+=(r**n/factn)
                    p00+=((r**C)/factC)*((C*mu)/C*(mu-lamda))
                p0=1/p00
                Lq0=((r**(C+1))/C)/(factC*((1-r/C)**2))
                Lq1=Lq0*p0
                Lq.set(round(Lq1,3))
                Wq.set(round((Lq1/lamda),3))
                w0=(Lq1/lamda)+(1/mu)
                W.set(round(w0,3))
                L.set(round((Lq1+r),3))
                
            frame2.grid_forget()
            frame4.grid_forget()
            frame5.grid_forget()
            frame3.grid(row = 4, column=1,columnspan=4,padx=160)
            frame3e1=Entry(frame3,width=30,borderwidth=4,bg="#202024",fg="white")
            frame3e1.grid(row=1,column=2,columnspan=1)
            frame3e2=Entry(frame3,width=30,borderwidth=4,bg="#202024",fg="white")
            frame3e2.grid(row=1,column=4,columnspan=1)
            frame3e3=Entry(frame3,width=30,borderwidth=4,bg="#202024",fg="white")
            frame3e3.grid(row=2,column=3,columnspan=1,pady=10)
            e1=Label(frame3,text="λ",bg="#202024",fg="white").grid(row=1,column=1)
            e2=Label(frame3,text="μ",bg="#202024",fg="white").grid(row=1,column=3,sticky=E)
            e3=Label(frame3,text="c",bg="#202024",fg="white").grid(row=2,column=2,sticky=E)
            btn0=Button(frame3,text="Calculate",padx=30,command=MMC)
            btn0.grid(row = 3, column=3,columnspan=1,pady=6)
            e4=Label(frame3,text="L",bg="#202024",fg="white").grid(row=4,column=1,sticky=E)
            e5=Label(frame3,text="Lq",bg="#202024",fg="white").grid(row=4,column=3,sticky=E)
            e6=Label(frame3,text="W",bg="#202024",fg="white").grid(row=5,column=1,sticky=E)
            e7=Label(frame3,text="Wq",bg="#202024",fg="white").grid(row=5,column=3,sticky=E)
            e8=Label(frame3,textvariable=L,padx=80,bg="#327ba8").grid(row=4,column=2,pady=5)
            e9=Label(frame3,textvariable=Lq,padx=80,bg="#327ba8").grid(row=4,column=4,pady=5)
            e10=Label(frame3,textvariable=W,padx=80,bg="#327ba8").grid(row=5,column=2,pady=5)
            e11=Label(frame3,textvariable=Wq,padx=80,bg="#327ba8").grid(row=5,column=4,pady=5)

        elif(clicked2.get()=="M/M/1/K") :
            L = IntVar()
            Lq= IntVar()
            W = IntVar()
            Wq= IntVar()
            def MM1K() :
                global lamda
                global mu
                global K
                lamda=float(frame4e1.get())
                mu=float(frame4e2.get())
                K=int(frame4e4.get())
                roh = lamda/mu
                global Pk
                global l0
                if(roh==1):
                    
                    L.set(round(K/2,3))
                    l0=K/2
                    Pk=1/(K+1)
                else :
                    l0=roh*(((1-((K+1)*(roh**K)))+(K*(roh**(K+1))))/((1-roh)*(1-(roh**(K+1)))))
                    L.set(round(l0,3))
                    Pk=(roh**K)*((1-roh)/(1-(roh**(K+1))))
                lamdadash=lamda*(1-Pk) 
                w0=l0/lamdadash
                W.set(round(w0,3))
                Wq.set(round(w0-1/mu,3))
                Lq.set(round(lamdadash*(w0-(1/mu)),3))

            frame2.grid_forget()
            frame3.grid_forget()
            frame5.grid_forget()
            frame4.grid(row = 4, column=1,columnspan=4,padx=160)
            frame4e1=Entry(frame4,width=30,borderwidth=4,bg="#202024",fg="white")
            frame4e1.grid(row=1,column=2,columnspan=1)
            frame4e2=Entry(frame4,width=30,borderwidth=4,bg="#202024",fg="white")
            frame4e2.grid(row=1,column=4,columnspan=1)
            frame4e3=Entry(frame4,width=30,borderwidth=4,state=DISABLED,disabledbackground="#403f4a")
            frame4e3.grid(row=2,column=2,columnspan=1,pady=10)
            frame4e4=Entry(frame4,width=30,borderwidth=4,bg="#202024",fg="white")
            frame4e4.grid(row=2,column=4,columnspan=1,pady=10)
            e1=Label(frame4,text="λ",bg="#202024",fg="white").grid(row=1,column=1)
            e2=Label(frame4,text="μ",bg="#202024",fg="white").grid(row=1,column=3,sticky=E)
            e3=Label(frame4,text="c",bg="#202024",fg="white").grid(row=2,column=1,sticky=E)
            e4=Label(frame4,text="K",bg="#202024",fg="white").grid(row=2,column=3,sticky=E)
            btn0=Button(frame4,text="Calculate",padx=30,command=MM1K)
            btn0.grid(row = 3, column=3,columnspan=1,pady=6)
            e4=Label(frame4,text="L",bg="#202024",fg="white").grid(row=4,column=1,sticky=E)
            e5=Label(frame4,text="Lq",bg="#202024",fg="white").grid(row=4,column=3,sticky=E)
            e6=Label(frame4,text="W",bg="#202024",fg="white").grid(row=5,column=1,sticky=E)
            e7=Label(frame4,text="Wq",bg="#202024",fg="white").grid(row=5,column=3,sticky=E)
            e8=Label(frame4,textvariable=L,padx=80,bg="#327ba8").grid(row=4,column=2,pady=5)
            e9=Label(frame4,textvariable=Lq,padx=80,bg="#327ba8").grid(row=4,column=4,pady=5)
            e10=Label(frame4,textvariable=W,padx=80,bg="#327ba8").grid(row=5,column=2,pady=5)
            e11=Label(frame4,textvariable=Wq,padx=80,bg="#327ba8").grid(row=5,column=4,pady=5)

        elif(clicked2.get()=="M/M/C/K") :
            L = IntVar()
            Lq= IntVar()
            W = IntVar()
            Wq= IntVar()
            def MMCK() :
                global lamda
                global mu
                global C
                global K
                lamda=float(frame5e1.get())
                mu=float(frame5e2.get())
                r = lamda/mu
                C=int(frame5e3.get())
                K=int(frame5e4.get())
                roh=r/C
                factC=1
                for fac in range(1,C+1):
                    factC=factC*fac
                factK=1
                for fac in range(1,K+1):
                    factK=factK*fac
                l00=0
                for B in range(1,C):
                    factb=1
                    for V in range(1,B+1):
                        factb=factb*V
                    l00+=(C-B)*(r**B)/factb
                p00=0
                if(roh!=1):
                    for n in range(0,C):
                        factn=1
                        for f in range(1,n+1):
                            factn=factn*f
                        p00+=(r**n/factn)
                    p00+=((r**C)/factC)*((1-(roh**(K-C+1)))/(1-roh))
                else :
                    for n in range(0,C):
                        factn=1
                        for f in range(1,n+1):
                            factn=factn*f
                        p00+=(r**n/factn)
                    p00+=((r**C)/factC)*(K-C+1)
                p0=1/p00
                Lq0=((roh*(r**C)*p0)/(factC*((1-roh)**2)))*(1-(roh**(K-C+1))-((1-roh)*(K-C+1)*roh**(K-C)))
                Pk=((lamda**K)/(((C**(K-C))*factC)*(mu**K))*p0)
                lamdadash=lamda*(1-Pk) 
                Lq.set(round(Lq0,3))
                L0=Lq0+C-(p0*(l00))
                L.set(round(L0,3))
                Wq.set(round((Lq0/lamdadash),3))
                W.set(round((L0/lamdadash),3))
                


            frame2.grid_forget()
            frame3.grid_forget()
            frame4.grid_forget()
            frame5.grid(row = 4, column=1,columnspan=4,padx=160)
            frame5e1=Entry(frame5,width=30,borderwidth=4,bg="#202024",fg="white")
            frame5e1.grid(row=1,column=2,columnspan=1)
            frame5e2=Entry(frame5,width=30,borderwidth=4,bg="#202024",fg="white")
            frame5e2.grid(row=1,column=4,columnspan=1)
            frame5e3=Entry(frame5,width=30,borderwidth=4,bg="#202024",fg="white")
            frame5e3.grid(row=2,column=2,columnspan=1,pady=10)
            frame5e4=Entry(frame5,width=30,borderwidth=4,bg="#202024",fg="white")
            frame5e4.grid(row=2,column=4,columnspan=1,pady=10)
            e1=Label(frame5,text="λ",bg="#202024",fg="white").grid(row=1,column=1)
            e2=Label(frame5,text="μ",bg="#202024",fg="white").grid(row=1,column=3,sticky=E)
            e3=Label(frame5,text="c",bg="#202024",fg="white").grid(row=2,column=1,sticky=E)
            e4=Label(frame5,text="K",bg="#202024",fg="white").grid(row=2,column=3,sticky=E)
            btn0=Button(frame5,text="Calculate",padx=30,command=MMCK)
            btn0.grid(row = 3, column=3,columnspan=1,pady=6)
            e4=Label(frame5,text="L",bg="#202024",fg="white").grid(row=4,column=1,sticky=E)
            e5=Label(frame5,text="Lq",bg="#202024",fg="white").grid(row=4,column=3,sticky=E)
            e6=Label(frame5,text="W",bg="#202024",fg="white").grid(row=5,column=1,sticky=E)
            e7=Label(frame5,text="Wq",bg="#202024",fg="white").grid(row=5,column=3,sticky=E)
            e8=Label(frame5,textvariable=L,padx=80,bg="#327ba8").grid(row=4,column=2,pady=5)
            e9=Label(frame5,textvariable=Lq,padx=80,bg="#327ba8").grid(row=4,column=4,pady=5)
            e10=Label(frame5,textvariable=W,padx=80,bg="#327ba8").grid(row=5,column=2,pady=5)
            e11=Label(frame5,textvariable=Wq,padx=80,bg="#327ba8").grid(row=5,column=4,pady=5)

    

        
confirm_btn=PhotoImage(file=Path(__file__).with_name('btn.png'))
img_label=Label(image=confirm_btn)
btn=Button(root,image=confirm_btn,command=frames,padx=60,borderwidth=0)
btn.grid(row=3,column=1, columnspan=5)

root.mainloop()