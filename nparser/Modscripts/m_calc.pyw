from PyQt5.QtWidgets import QMainWindow,QPushButton,QLineEdit,QGridLayout,QVBoxLayout,QApplication,QWidget,QMessageBox,QShortcut,QSizePolicy
from PyQt5 import QtWidgets,QtCore
from PyQt5.QtGui import QKeySequence
from PyQt5.QtCore import Qt
import sys
import itertools as it
import nparser.regparser as rp
import re
import json 
# TODO: Do kalkulačky používat format=True, jinak format=False (aby se s tím dalo regulérně pracovat)
# TODO: Mít nastavení jesli připouštět řešení z C\R
# TODO: Nejdřív zjistit,jaké jsou v stringu symbolu a *veškeré prohledávání* pak provádět jen přes tyto symboly!    

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

class Calc(QMainWindow):
    def __init__(self,parent=None):
        super().__init__(parent)
        self.init_default_settings()
        self._centralWidget=QWidget()
        self.setWindowTitle("PyCalc")
        self.setCentralWidget(self._centralWidget)
        self.mainLayout=QVBoxLayout()
        self.error_dialog = QtWidgets.QErrorMessage(self) 
        self.load_settings()
        #self.setFixedSize(280, 235)
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
        self.display=StretchedLineEdit(self.settings["design"]["display_font_height_ratio"])
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
        self.cols=6
        self.display.returnPressed.connect(lambda : self._writeDisplay("="))
        self.error_dialog.setWindowModality(QtCore.Qt.WindowModal)
        self._setupButtons()
        self.parser=self.init_parser()
        self.memory=[]
        self.mempos=-1
        self.undoShortcut=QShortcut(QKeySequence("Ctrl+M"),self)
        self.redoShortcut=QShortcut(QKeySequence("Ctrl+N"),self)
        self.clrMemShortcut=QShortcut(QKeySequence("Ctrl+,"),self)
        self.undoShortcut.activated.connect(lambda dir="backwards":self.display_memory(direction=dir))
        self.redoShortcut.activated.connect(lambda dir="forward":self.display_memory(direction=dir))
        self.clrMemShortcut.activated.connect(self.clear_memory)
        self.plotShortcut=QShortcut(QKeySequence("Ctrl+P"),self)
        self.mplotShortcut=QShortcut(QKeySequence("Ctrl+O"),self)
        self.plotShortcut.activated.connect(self.launch_plot)
        self.mplotShortcut.activated.connect(lambda: self.launch_plot(plot2d=True))

    def init_default_settings(self):
        self.default_settings={"design":{"dimensions":[300,235],"button_style": "background-color: None","display_style":"background-color: None",                                   "button_font_height_ratio": 0.3,"display_font_height_ratio": 0.4,"display_to_keyboard_percentage":18},
                                "parser":{"implicit_conversion":True,"implicit_ops":True,
                                        "fractions":True,"solve_w_complex":True,"solformat":True},
                                "draw_division":{"use_colors":True,"color_mode":1,"size":18},"draw_parse_tree":{"size0":20}}
    def _setupButtons(self):
        positions=it.product(range(self.rows),range(self.cols))
        labels=[["0",".","C","CC","=","Draw"],[str(i) for i in range(1,4)]+["(",")",":"],
                [str(i) for i in range(4,7)]+["+","-","x"],[str(i) for i in range(7,10)]+["*","/","^"],
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
            ##print(f"{row},{col}: {self.buttons[row,col].text()}")
        for name,button in self.buttons.items():
            #button.setFixedSize(40,40)
            button.setMinimumHeight(40)
            button.setMinimumWidth(40)
            try:
                but_style=self.settings["design"]["button_style"]
                if but_style is not None:
                    button.setStyleSheet(but_style)
            except KeyError:
                pass
            f= lambda textarg=button.text(): textarg
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
                        text=self.parser.sparse(ct)
                except  (ValueError,NotImplementedError,TypeError,ZeroDivisionError,AttributeError) as e:
                    self.error_dialog.showMessage(str(e))
                    #self.show_error_box(str(e))
                    text=ct
                cp=len(text)
                #text=str(eval(ct))
            else:
                text=""
        elif text=="Draw":
            ct=self.display.text()
            cp=self.display.cursorPosition()
            if ct:
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
                
                text=ct

                #text=str(eval(ct))
            else:
                pass
        elif text=="CC":
            text=""
        else:
            cp=self.display.cursorPosition()
            old_text=self.display.text()
            if text=="C":
                if cp==0:
                    text=old_text[1:]
                else:
                    text=old_text[:cp-1]+old_text[cp:]
                cp=cp
            else:
                text=old_text[:cp]+text+old_text[cp:]
                cp=cp+1
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
        if text not in self.memory and text.strip()!="":
            self.memory.append(text)
    def display_memory(self,direction="backwards"):
        try:
            text=self.memory[self.mempos]
            if direction=="backwards":
                if self.mempos==-len(self.memory):
                    self.mempos=-1
                else:
                    self.mempos-=1
            else:
                if self.mempos==len(self.memory)-1:
                    self.mempos=0
                else:
                    self.mempos+=1
            
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
        args=("from","to","style","pden","grid","polar","cont")
        llim=None;ulim=None;style="";pden=None
        from_res,to_res,style_res,pden_res,grid_res,polar_res,cont_res=[list(self.parser.find_all(text,pat)) for pat in args]
        if len(from_res)>1:
            raise ValueError("Please specify only a single lower bound (i.e. 'from' can only appear once)")
        if len(to_res)>1:
              raise ValueError("Please specify only a single upper bound (i.e. 'to' can only appear once)")
        if len(style_res)>1:
              raise ValueError("Please specify only a single style keyword")
        if len(grid_res)>1:
              raise ValueError("Please specify only a single grid keyword")
        if len(pden_res)>1:
              raise ValueError("Please specify only a single pden keyword")
        if len(polar_res)>1:
              raise ValueError("Please specify only a single polar keyword")
        if len(cont_res)>1:
              raise ValueError("Please specify only a single cont keyword")    

        #print("Got here")
        text=self.parser.strip_brackets(text,aggresive=False)
        #print(text)
        indices=[text.find(el) for el in args]
        found=any(index!=-1 for index in indices)
        
     
        
            
        use_grid=None
        polar=None
        contours=None
        
        #print(found)
        if found:
            minind=min(ind for ind in indices if ind!=-1)
            #print("Got here")
            main,bounds=text[:minind],text[minind:]
            ##print("#printing",main,bounds)
            #print(0)
            dlimit="(?:"+"|".join(re.escape(token) for token in args)+"|$)"
            #print(dlimit)
            #print("bounds",bounds)
            #from_part=re.findall("from(.*?)(?:to|from|style|pden|$)",bounds)
            #to_part=re.findall("to(.*?)(?:from|to|$)",bounds)
            from_part=re.findall(f"from(.*?){dlimit}",bounds)
            #from_part=re.findall("from(.*?)(?:to|from|style|pden|$)",bounds)
            #print("From_part",from_part)
            to_part=re.findall(f"to(.*?){dlimit}",bounds)
            #to_part=re.findall("to(.*?)(?:to|from|style|pden|$)",bounds)
            style_part=re.findall(f"style(.*?){dlimit}",bounds)
            #style_part=re.findall("style(.*?)(?:to|from|style|pden|$)",bounds)
            pden_part=re.findall(f"pden(.*?){dlimit}",bounds)
            #pden_part=re.findall("pden(.*?)(?:to|from|style|pden|$)",bounds)
            grid_part=re.findall(f"grid(.*?){dlimit}",bounds)
            polar_part=re.findall(f"polar(.*?){dlimit}",bounds)
            cont_part=re.findall(f"cont(.*?){dlimit}",bounds)
            if from_part:
                try:
                    llim=float(from_part[0])
                except (ValueError,TypeError):
                    try:
                        res=self.parser.parse(from_part[0],draw=False)
                        try:
                            llim=[float(el) for el in res]
                        except (ValueError,TypeError):
                            llim=float(res)
                    except (TypeError,ValueError):
                        raise ValueError(f"The lower bound '{str(res)}' can't be interpreted as a number!")
            if to_part:
                try:
                    ulim=float(to_part[0])
                except (ValueError,TypeError):
                    try:
                        res=self.parser.parse(to_part[0],draw=False)
                        try:
                            ulim=[float(el) for el in res]
                        except (ValueError,TypeError):
                            ulim=float(res)
                    except (TypeError,ValueError):
                        raise ValueError(f"The upper bound '{str(res)}' can't be interpreted as a number!")

            
            if style_part:
                style=style_part[0].strip()
            
            if pden_part:
                pden=int(pden_part[0])
            if grid_part:
                g=grid_part[0].strip()
                if g=="on" or g=="" or g.strip()=="1":
                    use_grid=True
                elif g=="off" or  g.strip()=="0":
                    use_grid=False
                else:
                    raise ValueError("Specify 'grid' with either 'on'/'off' or '1','0'" )
            if polar_part:
                g=polar_part[0].strip()
                if g=="on" or g=="" or g.strip()=="1":
                    polar=True
                elif g=="off" or  g.strip()=="0":
                    polar=False
                else:
                    raise ValueError("Specify 'polar'  with either 'on'/'off' or '1'/'0'" )
            if cont_part:
                g=cont_part[0].strip()
                if g=="on" or g=="" or g.strip()=="1":
                    contours=True
                elif g=="off" or  g.strip()=="0":
                    contours=False
                else:
                    raise ValueError("Specify 'cont'  with either 'on'/'off' or '1'/'0'" )
        else:
            main=text
        #print("parsing ",main)
        res=self.parser.parse(main)
        
        if plot2d:
            reslist=[]
            if isinstance(res,rp.Vector):
                if any (isinstance(comp,rp.Vector) for comp in res):
                    reslist=list(res)
                else:
                    reslist=(res,)
            else:
                reslist=(res,)
            self.parser.plot2(*reslist,llim=llim,ulim=ulim,style=style,pden=pden,use_grid=use_grid,polar=polar,contours=contours)
        else:
            if not type(res) ==rp.Vector:
                res=(res,)
            else:
                vars=res.get_vars()
                if len(vars)==1 and next(iter(vars))=="t" and res.dim==2:
                    if ulim==None:
                        ulim=2*3.14

                    self.parser.plot(res,llim=llim,ulim=ulim,style=style,pden=pden,use_grid=use_grid,polar=polar,contours=contours)
                    return
            self.parser.plot(*res,llim=llim,ulim=ulim,style=style,pden=pden,use_grid=use_grid,polar=polar,contours=contours)
            
if __name__=="__main__": 
    app=QApplication(sys.argv)
    menu=Calc()
    menu.show()
    sys.exit(app.exec())