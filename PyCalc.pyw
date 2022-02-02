from PyQt5.QtWidgets import QMainWindow,QPushButton,QLineEdit,QGridLayout,QVBoxLayout,QApplication,QWidget,QMessageBox,QShortcut,QSizePolicy
from PyQt5 import QtWidgets,QtCore
from PyQt5.QtGui import QKeySequence,QIcon
from PyQt5.QtCore import Qt,QEvent
import sys
import itertools as it
import nparser.regparser as rp
import numpy as np
import matplotlib.pyplot as plt
import re
import json 
# TODO: Do kalkulačky používat format=True, jinak format=False (aby se s tím dalo regulérně pracovat)
# TODO: Mít nastavení jesli připouštět řešení z C\R (done)
# TODO: α#1+2*3 se interpretuje jako α, mělo by se to interpretovat jako chybný vstup   
# TODO: 5Δxyz vyhodí chybu, prepraser by měl za 5 doplnit * (mělo by stačit jen rozšířit doplňování na unární operátory) 
# TODO: Nahrazování indexů v draw_parse_tree pořád někdy ostraňuje závorky, chce to tu verzi z plotů
# TODO: při appendování do paměti posunout pozici o +1 - done
# TODO: Ošetřit overflow exception, keyboard interrupt v event (try)
# TODO: Přidat více tlačítek, tlačítko pro nastavení plotu
# TODO: Načítání ze souboru
# TODO: Hezčí vykreslování parsovacího stromu (a je tam třeba opravit vstup u unárních operátorů
# TODO: Možná dokonce těch sliců budeme chtít víc najednou... číslo => všechny slicy s tím číslem, dvojice čísel => slice, vektor dvojice čísel => všechny takové dvojice
# TODO: slice tedy i pro vektorová pole...pokud je slice None, tak ukázat všechny, jinak zvolený slice
# TODO: Mít správně labely pro ode, i vzhledem k slicům
# TODO: Mít dx,dy pro diferenciální veličiny je nepraktické, protože se musí pořád zadávat |dx;dy -> co mít korespondenci dx <->d, dy <->e ? 
# TODO: dx a dy zatím vyřešeno přidáním pevných proměnných dx a dy...
# TODO: Celé to kreslení konstant je nějaké divné...například když se vloží (explicitní) vektor, tak se  mnohokrát překreslí ta samá věc
# TODO: Maps: komplexní hodnoty
# TODO: BVP: Přidat možnost a) definovat proměnné b) definovat parametry
# TODO: Přidat možnost specifikovat "title"
# TODO: x_1*7 -> xxxxxxx (opravit, asi by stačilo buď u variablu to vůbec nepřidávat jako attached argument nebo tam zakázat indexování)
# nice: [b(k-1)]sin(t)-bsin[t(k-1)],[b(k-1)]cos(t)+bcos[t(k-1)] |b=2;k=0.25,0.5,1.25,5 from -6pi to 6pi
# s[b(k-1)]sin(t)-bsin[t(k-1)],[b(k-1)]cos(t)+bcos[t(k-1)] |s=-1,1;b=2;k=0.25 from -6pi to 6pi 
# regparser: "V <- x" zruší celou vektorovou funkci, měla by být možnost to vrátit zas zpátky...
# regparser: ode a bvp se chovají divně, když je moc málo proměnných, například tohle dá chybu (špatný slice): ode([t,0],(0,0)) -> Mělo by být vždy n+1 proměnných, kde n je dimenze počátečního bodu (a ideálně y až druhá proměnná, pokud n>=2)
# regparser: bvp...nedá se nějak lépe vymyslet specifikace těch okrajů? Třeba přes čísla argumentů schovaných do vektoru: ((1,1,0),(1,2,-20))
# dilatace času: (1/gamma*t,k) |gamma-> (1/√[1-(0.025*k)^2]) ; k=0..39  from 0 to 10 anim point 
# šikmý vrh s odporem vzduchu: ode([H(y)dx,H(y)dy,-knorm((dx,dy))sign(dx),g-knorm((dx,dy))sign(dy)],[0,1,10,10]) |g=-10;k=0,0.05..1 t1 5 tden 1000 anim point realtime 
# Rovnice na výpočet BMR: 11*m-3*a+272*g+777 |m=80;a=30;g=1
# TODO: více řádků
# TODO: možnost ukládat výrazy na HD
class StretchedButton(QPushButton):
    def __init__(self,lab,resfactor):
        self.resfactor=resfactor
        super().__init__(lab)

    def resizeEvent(self, evt):
        font = self.font()
        font.setPixelSize(self.height() * self.resfactor)
        self.setFont(font)
class StretchedLineEdit(QLineEdit):
    def __init__(self,resfactor):
        self.resfactor=resfactor
        super().__init__()
    def resizeEvent(self, evt):
        font = self.font()
        font.setPixelSize(self.height() * self.resfactor)
        self.setFont(font)
# New code ↓
class myLineEdit(StretchedLineEdit):
    def __init__(self,resfactor,parser):
        super().__init__(resfactor)
        self.parser=parser
        self.cards={"diamond":"♦","club":"♣","heart":"♥","spade":"♠"}
    def event(self,event):
        try:
            if event.type()==QEvent.KeyPress and (event.key()==Qt.Key_Space):
                cp=self.cursorPosition()
                text=self.text()
                subtext,rest=text[:cp],text[cp:]
                ind=subtext.find("\\")
                if ind!=-1:
                    if subtext[ind+1:cp] in("nabla","del"):
                        newsymb="∇"
                        newtext=subtext[:ind]+newsymb+rest
                    elif subtext[ind+1:cp]=="dot":
                        newsymb="·"	
                        newtext=subtext[:ind]+newsymb+rest
                    elif subtext[ind+1:cp]=="times":
                        newsymb="×"
                        newtext=subtext[:ind]+newsymb+rest
                    elif subtext[ind+1:cp]=="div":
                        newsymb="∇·"+"["
                        newtext=subtext[:ind]+newsymb+"]"+rest
                        #newtext=subtext[:ind]+newsymb+rest
                    elif subtext[ind+1:cp]=="rot":
                        newsymb="∇×"+"["
                        newtext=subtext[:ind]+newsymb+"]"+rest
                        #newtext=subtext[:ind]+newsymb+rest
                    elif subtext[ind+1:cp]=="grad":
                        newsymb="∇"+"["
                        newtext=subtext[:ind]+newsymb+"]"+rest
                        #newtext=subtext[:ind]+newsymb+rest
                    elif subtext[ind+1:cp]=="lap":
                        newsymb="Δ"+"[" 
                        newtext=subtext[:ind]+newsymb+"]"+rest 
                        #newtext=subtext[:ind]+newsymb+rest
                    elif subtext[ind+1:cp]=="sqrt":
                        newsymb="√"+"["
                        newtext=subtext[:ind]+newsymb+"]"+rest
                        #newtext=subtext[:ind]+newsymb+rest
                    elif subtext[ind+1:ind+4]=="der":
                        newsymb="der["
                        newtext=subtext[:ind]+newsymb+","+subtext[ind+4:cp]+"]"+rest
                    elif subtext[ind+1:cp] in self.parser.funcsymbols: # parser.funcsymbols...
                        newsymb=subtext[ind+1:cp]+"["
                        newtext=subtext[:ind]+newsymb+"]"+rest
                    
                    elif subtext[ind+1:ind+3]=="dd":
                        #newsymb="der["
                        newsymb="("
                        #newtext=subtext[:ind]+newsymb+","+subtext[ind+3:cp]+"]"+rest
                        newtext=subtext[:ind]+newsymb+")_["+subtext[ind+3:cp]+"]"+rest
                    elif subtext[ind+1:cp] in self.cards:
                        newsymb=self.cards[subtext[ind+1:cp]]
                        newtext=subtext[:ind]+newsymb+rest
                    else:
                        newsymb=rp.PrettyVariable.pretty(subtext[ind+1:cp])
                        newtext=subtext[:ind]+newsymb+rest
                    self.setText(newtext)
                    self.setCursorPosition(ind+len(newsymb))
                    return True
        
            return super().event(event)
        except KeyboardInterrupt:
            pass

class Calc(QMainWindow):
    def __init__(self,parent=None):
        super().__init__(parent)
        self.init_default_settings()
        self._centralWidget=QWidget()
        self.setWindowTitle("PyCalc")
        self.setWindowIcon(QIcon("pi.ico"))
        self.setCentralWidget(self._centralWidget)
        self.mainLayout=QVBoxLayout()
        self.error_dialog = QtWidgets.QErrorMessage(self) 
        self.error_dialog.setWindowTitle("Error")
        self.load_settings()
        #self.setFixedSize(280, 23(1((5)
        #self.setFixedSize(235, 235)
        
        
        if "dimensions" in self.settings["design"]:
            
            dimensions=self.settings["design"]["dimensions"]
            ##print(dimensions)
            if isinstance(dimensions,list):
                ##print("resizing",dimensions)
                self.resize(*dimensions)
            else:
                self.resize(300,300)
        else:
            self.resize(300,300) 
                
        
        ##print(type(self.parser))
        self._centralWidget.setLayout(self.mainLayout)
        # New code ↓
        self.parser=self.init_parser()
        #self.display=StretchedLineEdit(self.settings["design"]["display_font_height_ratio"]) 
        self.display=myLineEdit(self.settings["design"]["display_font_height_ratio"],self.parser)
        # New code ↑
        #self.display=StretchedLineEdit(self.settings["design"]["display_font_height_ratio"])
        try:
            self.display.setStyleSheet(self.settings["design"]["display_style"])
        except KeyError:
            pass
        self._setDisplay()
        self.buttons_layout=QGridLayout()
        self.mainLayout.addWidget(self.display,35)
        self.display.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding) 
        self.mainLayout.addLayout(self.buttons_layout,int(35*(100/self.settings["design"]["display_to_keyboard_percentage"]-1)))
        self.rows=4
        self.cols=7
        self.display.returnPressed.connect(lambda : self._writeDisplay("="))
        self.error_dialog.setWindowModality(QtCore.Qt.WindowModal)
        self._setupButtons()
        self.memory=[]
        self.mempos=-1
        self.undoShortcut=QShortcut(QKeySequence("Ctrl+M"),self)
        self.redoShortcut=QShortcut(QKeySequence("Ctrl+N"),self)
        self.clrMemShortcut=QShortcut(QKeySequence("Ctrl+-"),self)
        self.undoShortcut.activated.connect(lambda dir="backwards":self.display_memory(direction=dir))
        self.redoShortcut.activated.connect(lambda dir="forward":self.display_memory(direction=dir))
        self.clrMemShortcut.activated.connect(self.clear_memory)
        self.plotShortcut=QShortcut(QKeySequence("Ctrl+P"),self)
        self.mplotShortcut=QShortcut(QKeySequence("Ctrl+I"),self)
        self.plotShortcut.activated.connect(self.launch_plot)
        self.mplotShortcut.activated.connect(lambda: self.launch_plot(plot2d=True))
        self.clstShortcut=QShortcut(QKeySequence("Ctrl+Q"),self)
        self.clstShortcut.activated.connect(lambda: plt.close("all"))
    def init_default_settings(self):
        self.default_settings={"design":{"dimensions":[300,235],"button_style": "background-color: None","display_style":"background-color: None",  "button_font_height_ratio": 0.3,"display_font_height_ratio": 0.4,"display_to_keyboard_percentage":18},
                                "parser":{"implicit_conversion":True,"implicit_ops":True,"restrict_symb": False,
                                        "fractions":True,"solve_w_complex":True,"solformat":True,"backend":"Qt5Agg"},
                                "draw_division":{"use_colors":True,"color_mode":1,"size":18},"draw_parse_tree":{"size0":20}}
    def _setupButtons(self):
        positions=it.product(range(self.rows),range(self.cols))
        labels=[["0",".","C","CC","=","Draw","2D plot"],[str(i) for i in range(1,4)]+["()"," | ",":","Plot"],
                [str(i) for i in range(4,7)]+["+","-","x","to"],[str(i) for i in range(7,10)]+["*","/","^","from"],
        ] 
        self.buttons={}
        for row,col in positions: 
            label=labels[row][col]
            button=self.buttons[row,col]=StretchedButton(label,self.settings["design"]["button_font_height_ratio"])
            button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding) 
            self.buttons_layout.addWidget(self.buttons[row,col],self.rows-row,col)
            if label=="CC":
                self.buttons[row,col].setShortcut(QKeySequence("Ctrl+Backspace"))
            elif label=="Draw":
                self.buttons[row,col].setShortcut(QKeySequence("Ctrl+d"))
            elif label=="from":
                self.buttons[row,col].setShortcut(QKeySequence("Ctrl+f"))
            elif label=="to":
                self.buttons[row,col].setShortcut(QKeySequence("Ctrl+t"))
            elif label=="()":
                self.buttons[row,col].setShortcut(QKeySequence("Ctrl+b"))
            elif label==" | ":
                self.buttons[row,col].setShortcut(QKeySequence("Ctrl+w"))
            
            ##print(f"{row},{col}: {self.buttons[row,col].text()}")
        for button in self.buttons.values():
            #button.setFixedSize(40,40)
            button.setMinimumHeight(40)
            button.setMinimumWidth(40)
            try:
                but_style=self.settings["design"]["button_style"]
                if but_style is not None:
                    button.setStyleSheet(but_style)
            except KeyError:
                pass
            #f= lambda textarg=button.text(): textarg
            ##print(name,f(),button.text())
            
            button.clicked.connect(lambda _,textarg=button.text(): self._writeDisplay(text=textarg))
    def _setDisplay(self):
        #self.display.setFixedHeight(35)
        self.display.setMinimumHeight(35)
        self.display.setAlignment(Qt.AlignRight)
        self.display.setReadOnly(False)
        self.display.setFocus() 
    def _writeDisplay(self,text):
        cp=0
        if text=="=":
            ct=self.display.text().replace("\\","\\")
          #  self.error_dialog.showMessage(ct)
            #return "2"
            if ct:
                text=ct.lstrip()
                ##print(ct)
                try:
                    self.update_memory(text)
                    if text.lower().startswith("(s)"):
                        self.proc_and_save(text,redraw=False,save=True) 
                        text=ct
                    elif text.lower().startswith("(m)"):
                        
                        self.mempos=-1
                        text=self.proc_and_save(text,redraw=False,save=False) 
                    elif text.lower().startswith("plot"):
                       self.handle_plot(text[4:]) 
                    elif text.lower().startswith("(p)"):
                       
                       self.handle_plot(text[3:]) 
                    elif text.lower().startswith("mplot"):
                       self.handle_plot(text[5:],plot2d=True) 
                    elif text.lower().startswith("(mp)"):
                       
                       self.handle_plot(text[4:],plot2d=True) 
                    else:
                       
                        self.mempos=-1
                        res=self.parser.parse(ct)
                        
                        if type(res) in (rp.Equation,rp.Inequality,rp.SystemOfEqs):
                            res=res.eval()
                        if type(res)==rp.ImplicitVector and all(type(comp) in (rp.Equation,rp.Inequality,rp.SystemOfEqs) for comp in res):
                            res=res.eval()
                        text=str(res)
                except  (ValueError,NotImplementedError,TypeError,ZeroDivisionError,AttributeError) as e:
                    self.error_dialog.showMessage(str(e))
                    #self.show_error_box(str(e))
                    text=ct
                cp=len(text)
                #text=str(eval(ct))
            else:
                text=""
        elif text in ("Draw","Plot","2D plot"):
            ct=self.display.text()
            cp=self.display.cursorPosition()
            if ct:
                if text=="Draw":
                    ##print(ct)
                    text=ct.lstrip()
                    try:
                        self.update_memory(text)
                        if text.lower().startswith("(s)"):
                            self.proc_and_save(text,redraw=True,save=True) 
                        elif text.lower().startswith("(m)"):
                            self.proc_and_save(text,redraw=True,save=False) 
                        elif text.lower().startswith("plot"):
                           self.handle_plot(text[4:]) 
                        elif text.lower().startswith("mplot"):
                           self.handle_plot(text[5:],plot2d=True) 
                        elif text.lower().startswith("(mp)"):
                           self.handle_plot(text[4:],plot2d=True) 

                        else:
                            lind=text.find("|")
                            if lind==-1:
                                lind=len(text)
                            if ":" in text and any( ind<lind for ind in self.parser.find_all(text,":")):
                                self.parser.longdiv(text,draw=True,**self.settings["draw_division"])
                            else:
                                self.parser.parse(text,draw=True,**self.settings["draw_parse_tree"])

                    except  (ValueError,NotImplementedError,TypeError,ZeroDivisionError,AttributeError) as e:
                        self.error_dialog.showMessage(str(e))
                        #self.show_error_box(str(e))

                    

                    #text=str(eval(ct))
                elif text=="Plot":
                    self.launch_plot()
                elif text=="2D plot":
                    self.launch_plot(plot2d=True)
            else:
                pass
            text=ct
        elif text=="CC":
            text=""
        elif text=="()":
            cp=self.display.cursorPosition()
            old_text=self.display.text()
            text="("+old_text+")"
            cp=len(text)
        else:
            cp=self.display.cursorPosition()
            old_text=self.display.text()
            symb=text
            if text=="C":
                if cp==0:
                    text=old_text[1:]
                else:
                    text=old_text[:cp-1]+old_text[cp:]
                cp=cp
            elif (word:=text) in ("from","to"):
                if word in old_text:
                    fpos=old_text.find(word)
                    if f"{word} " in old_text:    
                        text=old_text
                    else:
                        text=old_text[:fpos+len(word)]+" "
                    cp=fpos+len(word)+1
                else:
                    text=old_text+" "+f"{word} "
                    cp=len(text)
            else:
                text=old_text[:cp]+text+old_text[cp:]
                cp=cp+len(symb)
        self.display.setText(text)
        self.display.setCursorPosition(cp)
        self.display.setFocus()

    
    def show_error_box(self,e):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText(e)
        msg.setWindowTitle("Error")
        msg.show()

    def init_parser(self):
        #p=rp.bParser(implicit=True,verbose=False)
        p=rp.Parser(verbose=False,**self.settings["parser"])
        return p
    def load_settings(self):
        try:
            with open("../calc_settings.json") as setfile:
                
                self.settings=json.loads(setfile.read())
                #print("success")
        except PermissionError as e:
            self.error_dialog(str(e))
        except (FileNotFoundError) as e:  
            #self.error_dialog.showMessage(str(e))
            self.error_dialog.showMessage("Couldn't find 'calc_settings.json', reverting to default settings")
            
            self.settings=self.default_settings
            try:
                with open("../calc_settings.json","w") as f:
                    json.dump(self.settings,f,indent=2)
            except PermissionError as pe:
                self.error_dialog.showMessage(str(pe))
        except (ValueError) as e:
            self.error_dialog.showMessage(str(e))
            #print("Loading of settings failed, if you wish to use default settings, delete calc_settings.json")
            self.settings=self.default_settings
        self.settings={**self.default_settings,**self.settings}
    def update_memory(self,text):
        if text not in self.memory and text.strip()!="" :
            self.memory.append(text)
           # print(self.memory,self.mempos)
    def display_memory(self,direction="backwards"):
        try:
            if direction=="backwards":
                text=self.memory[self.mempos]
                if self.mempos==-len(self.memory):
                    self.mempos=-1
                else:
                    self.mempos-=1
            else:
                if self.mempos==-1:
                    self.mempos=-len(self.memory)
                else:
                    self.mempos+=1
                text=self.memory[self.mempos]
            #print(self.memory,self.mempos)
            if text==self.display.text() and len(self.memory)>1:
                self.display_memory(direction=direction)
            else:
                self.display.setText(text)
                self.display.setCursorPosition(len(text))
            
        except IndexError:
            self.display.setText("")
        
    def clear_memory(self): 
        self.memory=[]
        self.mempos=-1
    
    def proc_and_save(self,text,save=True,redraw=False):
        #with open("text.txt","w") as f:
        #    f.write("Ahoj")
        text=text[3:]
        if "|" in text:
            main_part,subs,*_=text.split("|")
            subs=" | "+subs 
        else:
            main_part,subs=text,""
        
        masked_text=self.parser.mask_text_im(main_part)
        def comma_seperated():
            last=-1
            for cind in self.parser.find_all(masked_text,","):
                yield last+1,cind
                last=cind
            yield last+1,None
        #commas=[-1]+list(self.parser.find_all(masked_text,","))
        #exprs=[main_part[i+1:j] for i,j in zip(commas,commas[1:]+[None]) ]
        exprs=[main_part[i:j] for i,j in comma_seperated() ]
        ##print(exprs)
        if save==redraw==False:

            reslist="(M) "+",".join([self.parser.sparse(expr+subs) for expr in exprs ])
            return reslist
        #    #print(exprs) 
        else:
            
            lind=(expr.find("|") for expr in exprs)
            norm_exprs=[]
            div_exprs=[]
            for expr in exprs:
                if ":" in expr:
                    lind=expr.find("|")
                    if lind==-1:
                        lind=len(text)
                    if  any( ind<lind for ind in self.parser.find_all(expr,":")):
                        div_exprs.append(expr)
                        continue
                norm_exprs.append(expr)
        #    print("1)",exprs,norm_exprs,div_exprs)
            if save: 

                try:
                    self.parser.draw_exprs(*(el+subs for el in norm_exprs),**self.settings["draw_parse_tree"])
                    self.parser.ldiv_and_save(*(el+subs for el in div_exprs),**self.settings["draw_division"])
                except (FileNotFoundError,PermissionError) as e:
                    self.error_dialog.showMessage(str(e))
                else:
                    QMessageBox.about(self,"Success","The figures were saved to the folder '../SavedFigures/'. ")
                    #self.error_dialog.showMessage("The figures were saved to the folder 'SavedFigures/'. ")
            if redraw:
                for expr in div_exprs:
                    self.parser.longdiv(expr+subs,draw=True,**self.settings["draw_division"])
                for expr in norm_exprs: 
                    self.parser.parse(expr+subs,draw=True,**self.settings["draw_parse_tree"])
    def launch_plot(self,plot2d=False):
        text=self.display.text().strip()
        self.update_memory(text)
        if text.lower().startswith("plot"):
            text=text[4:]
        if text.lower().startswith("mplot"):
            text=text[5:]
            plot2d=True
        elif text.lower().startswith("(p)"):
            text=text[3:]
        elif text.lower().startswith("(mp)"):
            text=text[4:]
            plot2d=True
        try:
            self.handle_plot(text,plot2d=plot2d)
        except  (ValueError,NotImplementedError,TypeError,ZeroDivisionError,AttributeError) as e:
                    self.error_dialog.showMessage(str(e))
    def handle_plot(self,text,plot2d=None):
        if plot2d is None:
            if text.startswith("2"):
                plot2d=True
                text=text[1:]
            else:
                plot2d=False
        args="from","to","style","pden","grid","polar","cont","hold","anim","equal","save","realtime","t0","t1","tden","slice","reim","legend","data","scale","loglog","hline","vline"
        llim=None;ulim=None;t0=None;t1=None;style="";pden=None;tden=None;use_grid=None;polar=None;contours=None;hold=None
        anim=None;point=None;equal=None;save=None;realtime=None;filename=None;v_slice=None;realimag=False;legend=True;get_data=False;scale=False;loglog=False;hlines=[];vlines=[]

        found={pat:list(self.parser.find_all(text,pat)) for pat in args}
        
        if len(found["from"])>1:
            raise ValueError("Please specify only a single lower bound (i.e. 'from' can only appear once)")
        if len(found["to"])>1:
              raise ValueError("Please specify only a single upper bound (i.e. 'to' can only appear once)")
        for word,res in found.items():
            if len(res)>1:
              raise ValueError(f"Please specify only a single '{word}' keyword")   
        
        #text,_=self.parser.strip_brackets(text,aggresive=False)
        indices=[text.find(el) for el in args]
        found=any(index!=-1 for index in indices)
        
        if found:
            minind=min(ind for ind in indices if ind!=-1)
            main,bounds=text[:minind],text[minind:]
            dlimit="(?:"+"|".join(re.escape(token) for token in args)+"|$)"
            parts={arg:re.findall(f"{arg}(.*?){dlimit}",bounds) for arg in args}
            """parts["from"]=re.findall(f"from(.*?){dlimit}",bounds)
            parts["to"]=re.findall(f"to(.*?){dlimit}",bounds)
            parts["style"]=re.findall(f"style(.*?){dlimit}",bounds)
            parts["pden"]=re.findall(f"pden(.*?){dlimit}",bounds)
            parts["grid"]=re.findall(f"grid(.*?){dlimit}",bounds)
            parts["polar"]=re.findall(f"polar(.*?){dlimit}",bounds)
            parts["cont"]=re.findall(f"cont(.*?){dlimit}",bounds)
            parts["hold"]=re.findall(f"hold(.*?){dlimit}",bounds)
            parts["anim"]=re.findall(f"anim(.*?){dlimit}",bounds)
            parts["equal"]=re.findall(f"equal(.*?){dlimit}",bounds)
            parts["save"]=re.findall(f"save(.*?){dlimit}",bounds)
            """
            if parts["from"]:
                try:
                    llim=float(parts["from"][0])
                except (ValueError,TypeError):
                    try:
                        res=self.parser.parse(parts["from"][0],draw=False)
                        try: 
                            llim=[float(el) for el in res]
                        except (ValueError,TypeError):
                            llim=float(res)
                    except (TypeError,ValueError):
                        raise ValueError(f"The lower bound '{str(res)}' can't be interpreted as a number!")
            if parts["to"]:
                try:
                    ulim=float(parts["to"][0])
                except (ValueError,TypeError):
                    try:
                        res=self.parser.parse(parts["to"][0],draw=False)
                        try:
                            ulim=[float(el) for el in res]
                        except (ValueError,TypeError):
                            ulim=float(res)
                    except (TypeError,ValueError):
                        raise ValueError(f"The upper bound '{str(res)}' can't be interpreted as a number!")
            if parts["t0"]:
                anim=True
                try:
                    t0=float(parts["t0"][0])
                except (ValueError,TypeError):
                    try:
                        expr=parts["t0"][0].strip()
                        if expr.startswith("="):
                            expr=expr[1:] 
                        res=self.parser.parse(expr,draw=False)
                        t0=float(res)
                    except (TypeError,ValueError):
                        raise ValueError(f"The lower bound '{str(res)}' can't be interpreted as a number!")
            if parts["t1"]:
                anim=True
                try:
                    t1=float(parts["t1"][0])
                except (ValueError,TypeError):
                    try:
                        expr=parts["t1"][0].strip()
                        if expr.startswith("="):
                            expr=expr[1:] 
                        res=self.parser.parse(expr,draw=False)
                        t1=float(res)
                    except (TypeError,ValueError):
                        raise ValueError(f"The upper bound '{str(res)}' can't be interpreted as a number!")
            if parts["tden"]: 
                anim=True
                try:
                    tden=int(parts["tden"][0])
                except (ValueError,TypeError):
                    try:
                        expr=parts["tden"][0].strip()
                        if expr.startswith("="):
                            expr=expr[1:] 
                        res=self.parser.parse(expr,draw=False)
                        tden=int(res)
                    except (TypeError,ValueError):
                        raise ValueError(f"The lower bound '{str(res)}' can't be interpreted as a number!")
            if parts["hline"]:
                try:
                    hlines=[float(parts["hline"][0])]
                except (ValueError,TypeError):
                    try:
                        res=self.parser.parse(parts["hline"][0],draw=False)
                        try: 
                            hlines=[float(el) for el in res]
                        except (ValueError,TypeError):
                            hlines=[float(res)]
                    except (TypeError,ValueError):
                        raise ValueError(f"The expression '{str(res)}' can't be interpreted as a number!")
            if parts["vline"]:
                try:
                    vlines=[float(parts["vline"][0])]
                except (ValueError,TypeError):
                    try:
                        res=self.parser.parse(parts["vline"][0],draw=False)
                        try: 
                            vlines=[float(el) for el in res]
                        except (ValueError,TypeError):
                            vlines=[float(res)]
                    except (TypeError,ValueError):
                        raise ValueError(f"The expression '{str(res)}' can't be interpreted as a number!")
            if parts["style"]:
                style=parts["style"][0].strip()
            if parts["pden"]:
                pden=int(parts["pden"][0])
            if parts["grid"]:
                g=parts["grid"][0].strip()
                if g=="on" or g=="" or g.strip()=="1":
                    use_grid=True
                elif g=="off" or  g.strip()=="0":
                    use_grid=False
                else:
                    raise ValueError("Specify 'grid' with either 'on'/'off' or '1','0'" )
            if parts["polar"]:
                g=parts["polar"][0].strip()
                if g=="on" or g=="" or g.strip()=="1":
                    polar=True
                elif g=="off" or  g.strip()=="0":
                    polar=False
                else:
                    raise ValueError("Specify 'polar'  with either 'on'/'off' or '1'/'0'" )
            if parts["reim"]:
                g=parts["reim"][0].strip()
                if g=="on" or g=="" or g.strip()=="1":
                    realimag=True
                elif g=="off" or  g.strip()=="0":
                    realimag=False
                else:
                    raise ValueError("Specify 'realimag'  with either 'on'/'off' or '1'/'0'" )
            if parts["data"]:
                get_data=True
            if parts["cont"]:
                g=parts["cont"][0].strip()
                if g=="on" or g=="":
                    contours=20
                
                elif g=="off":
                    contours=False
                else:
                    try:
                        contours=int(g)
                        if contours<0:
                            raise ValueError
                    except ValueError:
                        raise ValueError("Specify 'cont'  with either 'on'/'off' or a natural number " )
            if parts["hold"]:
                g=parts["hold"][0].strip()
                if g=="on" or g=="" or g.strip()=="1":
                    hold=True
                elif g=="off" or  g.strip()=="0":
                    hold=False
                else:
                    raise ValueError("Specify 'cont'  with either 'on'/'off' or '1'/'0'" )
            if parts["anim"]:                
                g=parts["anim"][0].strip()
               # print(g)
             #   print(g,g=="off",anim)
                if g=="on" or g=="" or g=="1":
                    anim=True
                elif g=="point":
                    anim=True
                    point=True
                elif g=="off" or  g=="0":
                    anim=False
                
                else:
                    raise ValueError("Specify 'cont'  with either 'on'/'off' or '1'/'0'" )
            if parts["realtime"]:
                g=parts["realtime"][0].strip()
                if g=="on" or g=="" or g.strip()=="1":
                    realtime=True
                    anim=True 
                elif g=="point":
                    realtime=True
                    anim=True
                    point=True
                elif g=="off" or  g.strip()=="0":
                    realtime=False
                
                else:
                    raise ValueError("Specify 'realtime'  with either 'on'/'off' or '1'/'0'" )
            if parts["equal"]:
                g=parts["equal"][0].strip()
             #   print(g,g=="off",anim)
                if g=="on" or g=="" or g=="1":
                    equal=True
                elif g=="off" or  g=="0":
                    equal=False
                else:
                    raise ValueError("Specify 'equal'  with either 'on'/'off' or '1'/'0'" )
            if parts["legend"]:
                g=parts["legend"][0].strip()
             #   print(g,g=="off",anim)
                if g=="on" or g=="" or g=="1":
                    legend=True
                elif g=="off" or  g=="0":
                    legend=False
                else:
                    raise ValueError("Specify 'legend'  with either 'on'/'off' or '1'/'0'" )
            if parts["slice"]:
                if parts["slice"][0].strip()=="all":
                    v_slice="all"
                else:
                    try:
                        v_slice=[int(s) for s in parts["slice"][0].strip().split(",")]
                    except ValueError:
                        raise ValueError("The input after slice must either be 'all' or  two integers seperated by a comma")
            if parts["scale"]:
                g=parts["scale"][0].strip()
             #   print(g,g=="off",anim)
                if g=="on" or g=="" or g=="1":
                    scale=True
                elif g=="off" or  g=="0":
                    scale=False
                else:
                    raise ValueError("Specify 'scale'  with either 'on'/'off' or '1'/'0'" )
            if parts["loglog"]:
                g=parts["loglog"][0].strip()
             #   print(g,g=="off",anim)
                if g=="on" or g=="" or g=="1":
                    loglog=True
                elif g=="off" or  g=="0":
                    loglog=False
                else:
                    raise ValueError("Specify 'loglog'  with either 'on'/'off' or '1'/'0'" )
            

            if parts["save"]:
                save=True
                g=parts["save"][0].strip()
             #   print(g,g=="off",anim)
                if g:
                    filename=g
                else:
                    filename=None
        else:
            main=text
        #print("parsing ",main)
        res=self.parser.parse(main)
        if get_data:
            self.parser.plotdata=[]
        if isinstance(res,(rp.Equation,rp.SystemOfEqs)):
            pass#res=res.eval()
        try:
            if plot2d:
                reslist=[]
                if isinstance(res,rp.ImplicitVector):
                    #if any (isinstance(comp,rp.Vector) for comp in res):
                    #    reslist=list(res)
                    #else:
                    #    reslist=(res,)
                    reslist=list(res)
                else:
                    reslist=(res,)
                self.parser.plot2(*reslist,llim=llim,ulim=ulim,style=style,pden=pden,use_grid=use_grid,polar=polar,contours=contours,equal=equal,save=save,filename=filename,
                realtime=realtime,anim=anim,t0=t0,t1=t1,tden=tden,realimag=realimag,v_slice=v_slice,legend=legend)
            else:
                try:
                    vars=res.get_vars()
                except AttributeError:
                    vars="x"
                if  len(vars)==1 and next(iter(vars))=="t":# and res.dim==2:
                        if llim is None:
                            llim=0
                        if ulim is None:
                            ulim=2*3.14
                if not isinstance(res,rp.ImplicitVector):
                    res=(res,)
                    
                        #if not any(isinstance(comp,rp.Vector)  for comp in res):
                        #    self.parser.plot(res,llim=llim,ulim=ulim,style=style,pden=pden,t0=t0,t1=t1,tden=tden,use_grid=use_grid,v_slice=slice,
                        #    polar=polar,contours=contours,hold=hold,anim=anim,pointanim=point,equal=equal,save=save,filename=filename,realtime=realtime)
                        #    return
                #print("res:",res)
                self.parser.plot(*res,llim=llim,ulim=ulim,style=style,pden=pden,t0=t0,t1=t1,tden=tden,use_grid=use_grid,polar=polar,contours=contours,
                                hold=hold,anim=anim,pointanim=point,equal=equal,save=save,filename=filename,realtime=realtime,v_slice=v_slice,realimag=realimag,legend=legend,scale=scale,loglog=loglog,hlines=hlines,vlines=vlines)
        except TypeError as e:
            print(e)
            self.error_dialog.showMessage("Error while plotting the graph (maybe unexpected complex values?)")
        except IndexError as e:
            print(e)
            self.error_dialog.showMessage("The indices of the slice are out of bounds!")
        except OverflowError as e:
            self.error_dialog.showMessage("Overflow error (the results are too large)! ")
        except FileNotFoundError as e:
            self.error_dialog.showMessage(str(e))
        if get_data and self.parser.plotdata:
            #raise ValueError(np.array(self.parser.plotdata))
            with open("plotdata.txt","w") as file:
                for ind,ar in enumerate(np.array(self.parser.plotdata)):
                    np.savetxt(file,ar.T)
                    if ind+1!=len(ar):
                        file.write("\n\n")  
    def closeEvent(self,event):
        plt.close("all")  
        super().closeEvent(event)

if __name__=="__main__": 
    #QApplication.setAttribute(Qt.AA_EnableHighDpiScaling)
    app=QApplication(sys.argv)
    menu=Calc()
    menu.show()
    sys.exit(app.exec())
    