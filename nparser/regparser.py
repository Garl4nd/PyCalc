from typing import Callable
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.figure import Figure
import math
import re
import time
from operator import mul
import functools as ft
import itertools as it
import numpy as np
from . import prime_factor
from . import solver
from . import QRLU
from . import gen_obj
from . import difeq as de
from .gen_obj import (SPolynomial,Polynomial,Monom,Vector,ImplicitVector,
Quotient,Variable,PrettyVariable,Complex,
Equation,Inequality,SystemOfEqs,Solution,GeneralObject)
from . import appr_polsolve
# TODO: Přidat do animací volbu fps - teď je možné jen rychlost 1:1 nebo jak vyjde podle velikosti časového kroků ((t1-t0)/tden)
# TODO: implicitni doplnění pro závorky a unární operátory? : (x-1)#2 -> (x-1)*#2 - done
# TODO: plotování komplexnich cisel do komplexní roviny - done
# TODO: x**5 hlásí dělení nulou, co za to může? - solved
# TODO: iterplot, bifurcation diagram - done
# TODO: zadání Lagrangiánu -> vyřešení pohybových rovnic a animace <3 
# TODO: -> by  konzistence mělo taky vracet jen jeden plochý vektor -> done
# TODO: "V(sin(xt),cos(xt)) anim" v plot(...) vyhodí chybu, přitom by to mohlo být vektorové pole. Pro 2D ploty to funguje jak má.
# TODO: Umožnit kreslení všech páru proměnných u řešení diferenciálních rovnic (všechny slicy), to znamená, že dif_eq_plot asi bude muset převzít zodpovědnost za vytváření, ukazování, ukládání a zavírání obrázků
# TODO: matmul,dot,cross by měly být obecné funkce, aby byly schopny přijímat obecné argumenty (vyhodnocování přes =). real a imag by měly mít obecné alternativy. 
# TODO: plot bodů by mohl umět ukazovat reálnou a imag část vs pořadí - done
# TODO: maticové výpočty - eigenvalues, inverses...nakonec už to implementovaný v solveru soustav lin. rovnic. 
# TODO: V(IV(x)) by mělo dát V(x) - done
# TODO: Když do plot2 vejde vektor konstantních vektorů, mají se nakreslit v rovině jako body? - ne, od toho je plot
# TODO: difeq by měla rovnou zahrnovat parametry, které se hned dosadí...tím by se to mohlo hodně zrychlit i zhezčit výpis
#TODO # mocnění na vektor byx mělo vrátit vektor by mělo vrátit vekttor mocnin - done
#TODO: # if not newsymb in text: provizorní? Ale if it ain't broke don't fix it.
# TODO: (x-t)^2+y^2=1 -> bylo by hezké dostat animaci posouvající se kružnice
# TODO: přidat možnost explicitně specifikovat, v jaké proměnné se má dělat animace. Návrh: "sin(xt) anim x" - done, ale místo anim x se píše avar x.
# TODO: H(t-u)real(sqrt((t-u)^2 - (x-u)^2+0i)) |u=1,1.2..10 realtime t0=0 t1=10   from 0 to 10  #ukazuje se jen druhá polovina?? Ale když se to nejdřív spočítá a až pak nakreslí, tak to jde
# TODO: Hezčí bifurkační funkce a lepší výstup (nebo vstup), aby to neukazovalo milion "BIfurc map with init point..." etc.
# TODO: Na netu nefunguje tden=500
# TODO: "?" hlásí chybu v re (tj. spadne Pycalc, protože to tam neošetřuju). 
# regparser: "V <- x" zruší celou vektorovou funkci, měla by být možnost to vrátit zas zpátky, nebo zakázat přiřazování k předdefinovaným proměnným

class Op:
    def __init__(self,symbol:str,pri: int,arity: int,action: Callable,
                identity:set = None,direction="right",acts_on="left",unary_overload=None):
        self.pri=pri
        self.arity=arity
        self.action=action
        self.symbol=symbol
        self.identity=identity
        self.direction=direction
        self.acts_on=acts_on
        self.unary_overload=unary_overload
        self._textual=False
        if self.identity:
            self.id_rep=next(iter(self.identity))
    
    def intr(self,*args):    
        
        if len(args)!=self.arity and self.arity!="*":
            raise ValueError("Spatny pocet argumentu!")
        else:
            ##print(self.action,args)
            return self.action(*args)
    def __repr__(self):
        return "Priority: {0}, arity: {1},interpretation: {2},direction: {3}, acts to the: {4}".format(
            self.pri,self.arity,self.action,self.direction,self.acts_on
        )
class Parser:


    def __init__(self,preload=True,preparse=True,fill_id=True,force_func_brackets=False,implicit_conversion=True,
                restrict_symb=True,verbose=True,implicit_ops=True,
                fractions=True,solve_w_complex=True,solformat=False,backend="TkAgg"):
        self.operators={}
        self.opsymbols=set()
        self.funcs={}
        self.funcsymbols=set()
        self.symbols=set([":=","->","|","where"])
        self.brackets=[]
        self.modifier="@"
        self.forbidden_symbols={"-","+"} # Plus,- a . pred cislem python chape jako cislo. ALternativa je mist maskovani kontrolovat, jestli se tyto symboly vyskytuji na zacatku a neuznavat je jako cisla.
        self.modified_symbols=set()
        self.implicit_op="*"
        self.implicit_conversion=implicit_conversion # prevadet automaticky nedefinovane symboly na polynomy?
        self.orig_text=None
        self.verbose=verbose
        self.preparse=preparse
        self.implicit_ops=implicit_ops
        self.fill_id=fill_id
        self.force_func_brackets=force_func_brackets
        self.restrict_symb=restrict_symb
        self.fractions=fractions
        self.solve_w_complex=solve_w_complex
        self.backend=backend
        self.plot_ind=0
        self.solformat=solformat
        self.max_symbol_len=0
        # emergency error handling 
        self.max_error_depth=10
        self.__cur_error_depth=0
        self.__error_text=None
        self.pollist=[]
        self.pops={}
        self.pfuncs={}
        self.def_func=""
        self.def_func2=""
        self.settings={"preparse":preparse,"verbose":verbose,"restrict_symb":restrict_symb,"fill_id":fill_id,"force_func_brackets":force_func_brackets,"implicit_conversion":implicit_conversion,
        "implicit_ops":implicit_ops,"fractions":fractions}
        if preload: 
            self.load_standard()
        if self.verbose:
            print(self)
        self.ac=self.AnimControl()
        self.animlist=[]
        self.plotdata=[]
    def show_settings(self,help=False):
       
        print("The parser is using the following settings: \n")
        if help:
            help_dict={"preparse": "Make minor changes to the string for easier processing (e.g. convert -- to +, convert *-4 to *(-4), check for correct placement of bracekts, etc.)",
                        "verbose": "Print additional information about the parsing process","restrict_symb": "Restrict implicit variable names to letters of the alphabet(s)",
                        "fill_id": "In case of missing arguments of binary operators, fill in their unit elements (eg. replace '*5+' with '1*5+0').",
                        "force_func_brackets":"Require that functional arguments are enclosed with brackets (e.g. 'sin 4 5' fails if this options is True, otherwise returns 5*sin(4)).",
                        "implicit_conversion":"If an undefined symbol is encountered, automatically convert all it's letters to variables",
                        "implicit_ops": "Automatically fill in a predifined operator between symbols and brackets (e.g. '4x' -> '4*x', '(4)(5)'->'(4)*(5)' " ,
                        "fractions" : "When solving equations, output the result as fractions (e.g. for '4x=1' return Q(1,4)=1/4 if True and 0.25 if False)"}
            for key,value in help_dict.items():
                print(f"{key} ({self.settings[key]}): {value}\n")
            print("\nTo temporarily change the settings, you can set the corresponding keyword argument when calling parse(): \n\np=Parser();\np.parse('4*5',verbose=False)"  
                    "\n\nTo change the settings permanently, either set the corresponding keyword arguments when initializing the parser:\n\np=Parser(verbose=False,force_func_brackets=True)"
                    "\n\nor set them manually: \n\np=Parser()\np.verbose=False")

        else:
            for key,value in self.settings.items(): 
                print(key,":",value)
            print("\nFor details on the settings, call this function again with help=True")
        
    def __str__(self):
        return "\nParser with operators {0}\n {1:>20} {2} \n {3:>20} {4} \n\n{5:>20}".format(
            [self.unmake_unique_text(el) for el in list(self.operators.keys())],"and functions",list(self.funcs.keys()),"and brackets",self.brackets,
            "Call 'show_settings()' to see settings.")
        
    def add_operator(self,op: str,pri: int,arity: int,action: Callable,identity:set =None,
                    direction="right",acts_on="left",implicit_conversion=False,unary_func=None,overwrite=False):
        if op in self.forbidden_symbols:
            op=self.make_unique(op)
        #print(type(self.operators),op,op in self.operators)
        if op in self.symbols:
            if overwrite:
                if op in self.operators:
                    self.delete_operator(op)
                elif op in self.funcs:
                    self.delete_function(op)
            else:
                raise ValueError("Operator uz je zaregistrovany! Pridejte 'overwrite=True' ")
        if arity==1:
            if acts_on=="right":
                direction="right"
            else:
                direction="left"
        if implicit_conversion:
            self.implicit_op=op
        if unary_func==None:
            uop=None
        else:
            uop=Op("",pri,1,unary_func[0],acts_on=unary_func[1])
        #print(uop)
        self.operators[op]=Op(op,pri,arity,action,identity=identity,direction=direction,acts_on=acts_on,unary_overload=uop)
        self.opsymbols.add(op)
        self.add_symbol(op)
        
        self._sort_ops()
        
    def add_symbol(self,symbol):
        self.symbols.add(symbol)
        self.max_symbol_len=max(map(len,self.symbols))

    def delete_operator(self,op):
        if op in self.forbidden_symbols:
            op=self.make_unique(op)
        try:
            del self.operators[op]
            self.opsymbols.remove(op)
        except KeyError:
            print(f"The parser does not have operator {op}! ")

    def delete_function(self,name):
        try:
            del self.funcs[name]
            self.funcsymbols.remove(name)
            self.symbols.remove(name)
        except KeyError:
            print(f"The parser does not have function {name}! ")
    def add_function(self,name: str,arity: int,action: Callable,identity:set =None,overwrite=False,textual=False):
        if name in self.symbols:
            if overwrite:
                if name in self.operators:
                    self.delete_operator(name)
                elif name in self.funcs:
                    self.delete_function(name)
            else:
                raise ValueError("Operator uz je zaregistrovany! Pridejte 'overwrite=True' ")
        self.funcs[name]=Op(name,10000,arity,action,identity)
        self.funcs[name]._textual=textual
        self.funcsymbols.add(name)
        self.add_symbol(name)
        self._sort_functions()
    def make_unique(self,symbol):
        nt=self.modifier+symbol+self.modifier
        self.modified_symbols.add((symbol,nt))
        return nt
    def unmake_unique_text(self,text):
        for el in self.modified_symbols:
            text=text.replace(el[1],el[0])
        return text
    def add_brackets(self,brackets: list):
        if len(brackets)!=2:
            raise ValueError("Must specify a pair of brackets!")
        self.brackets.append(brackets)
        self.add_symbol(brackets[0])
        self.add_symbol(brackets[1])

    def delete_bracekts(self,brackets: list):
        try:
            self.brackets.remove(brackets)
            self.symbols.remove(brackets[0])
            self.symbols.remove(brackets[1])
        except ValueError:
            print(f"The parser does not have brackets {brackets}")
            
    def _sort_op_key(self,op):
        if op[1].arity==0:
            return 10**7-len(op[0])
        else:
            return op[1].pri
            
    def _sort_ops(self):                      
        self.operators=dict(sorted(self.operators.items(),key=self._sort_op_key))
    def _sort_functions(self):
        self.funcs=dict(sorted(self.funcs.items(),key=lambda item:-len(item[0])))
    
    def find_all(self,text,p,direction="right",hungry=False):
        '''Yields all the positions of
        the pattern p in the string s.'''

        if direction=="right":
            if p=="":
                return
            i = text.find(p)
            while i != -1:
                yield i
                if hungry:
                    i = text.find(p, i+1)
                else:
                    i = text.find(p, i+len(p))
        else:
            if p=="":
                return
            i = text.rfind(p)
            while i != -1:
                yield i
                if hungry:
                    i = text.rfind(p, i+1)
                else:
                    i = text.rfind(p, i+len(p))


    def _filter_brackets(self,text:str):
        masked_inds=[[] for _ in self.brackets]
        #print(text)
        for bind,bracks in enumerate(self.brackets):
            left,right=bracks
            linds=[(pos,"l") for pos in self.find_all(text,left)]
            rinds=[(pos,"r") for pos in self.find_all(text,right)]
            inds=sorted(linds+rinds)
            lc=0
            rc=0
            for ind in inds:
                if lc==0:
                    lp=ind[0]
                if ind[1]=="l":
                    lc+=1
                else:
                     rc+=1
                if lc==rc:
                    rp=ind[0]
                    lc=0
                    rc=0
                    masked_inds[bind].append((lp,rp))
                
        return masked_inds
    
    def mask_text(self,text,brack_inds):
        masked_text=str(text)
        for ind_list in brack_inds:
            for lind,rind in ind_list:
                masked_text=masked_text[:lind]+" "*(rind-lind+1)+masked_text[rind+1:]
        #print("text: {0}, masked: {1}, brack_inds: {2}".format(text,masked_text,brack_inds))
        return(masked_text)


    def mask_text_im(self,text): #C onveniece routine, combines bracks and mask
        brack_inds=self._filter_brackets(text)
        return self.mask_text(text,brack_inds)

    def is_number(self,text):
        try:
            float(text)
        except ValueError:
            
            return False
        else:
            return True
    
    def convert_to_num(self,text):
        try:
            res=int(text)
        except ValueError:
            res=float(text)
        return res

    def strip_brackets(self,text,aggresive=True):
        
        #print("Stripping brackets from :",text)
        text=text.strip()
        found=False
        #print("Stripped spaces :",text)
        for left,right in self.brackets:
            if text.startswith(left) and text.endswith(right):
                
         #       print("Processing ",left,right)
          #      print(f"{text} starts with ",left)
                left_inds=list(self.find_all(text,left))
                right_inds=list(self.find_all(text,right))

                if (len(left_inds)!=len(right_inds)):
                    raise ValueError("Spatne uzavrene zavorky!")
                else:
                    #nums=(len(k for k in left_inds if k<i) for i in left_inds)
                    nums=(len([k for k in left_inds if k<ri])-(i+1) for i,ri in enumerate(right_inds[:-1])) #Pro kazdou pravou zavorku vrati kolik nalevo od ni prebyva levych zavorek. 
                    if any(num==0 for num in nums):
                        continue
                    else:
                        found=True        
                        if aggresive:
                            text,_=self.strip_brackets(text[len(left):-len(right)],aggresive=True)
                        else:
                            text=text[len(left):-len(right)]
                            break
                #print(f", stripping: {text}")
        if text.strip()=="":
            return "0",found
        else:
            return text,found

    
    def strip2(self,text,brack_inds):
        #print(f"strip text {text}",brack_inds)
        for ind,bind in enumerate(brack_inds):
            #print(ind,bind)0f
            if len(bind)==1:
                
                lind,rind=bind[0]
                len_left=len(self.brackets[ind][0])
                len_right=len(self.brackets[ind][1])
             #   print(len_left,len_right,lind,rind,text[len_left:rind])
                if lind==0 and rind==len(text)-len_right:
                    return text[len_left:rind]
        return text

    def check_forward_conflicts(self,text,symbol):
        for it_symbol in self.symbols:
            if len(it_symbol)>len(symbol):
                    if text.startswith(it_symbol):
                        return True
        return False


    def is_substring(self,text,my_pos,sub_len=1,let_pass={}):
        #print("checking",text,"|",text[my_pos:my_pos+sub_len],"...",my_pos)
        for symbol in self.symbols:
            if symbol in let_pass:
                continue
            for pos in self.find_all(symbol,text[my_pos:my_pos+sub_len]): # Najde pozici znaku v symbolu
                if my_pos-pos<0: # pozice znaku v symbolu je vetsi nez nase misto v textu (tudiz to nemuze byt match)
                    continue
                else:
                    if text[my_pos-pos:my_pos+len(symbol)-pos]==symbol:
                        return True,my_pos+len(symbol)-pos
        return False,None

    def contains_illegals(self,text,temp_ops=set()):
            letind=0
            #print("Checking",text)
            #input()
            while letind<len(text):
                #print(letind)
                letter=text[letind]
                if letter not in self.symbols.union({"."," ",",","|",}) and not letter.isdigit(): #is the letter missing from symbols?
                    found,new_ind=self.is_substring(text,letind)
                    if not found:
                        if self.implicit_conversion:
                            if letter.isalpha() or not self.restrict_symb:
                                #self.add_monvar(letter)
                                self.add_polvar(letter,temp_ops=temp_ops)
                                if self.verbose:
                                    print(f"Converting {letter} to variable")
                            else:
                                
                                raise ValueError(f"Can't implicitly convert the non-alphabetic symbol {letter} to variable. "
                            "To override this, set restrict_symb to 'False'")
                        else:
                            print("Found illegal letter ",letter,"which is not in",self.symbols)
                            return True,letter
                    else:
                        letind=new_ind
                else:
                    letind+=1
            return False,None

    def replace_unknown(self,replaced_symbol,new_symbol,text):
        letind=0
        #print("Checking",text)
        #input()
        end=len(text)
        newtext=""
        lo=len(replaced_symbol)
        ln=len(new_symbol)
        
       #print("n?","n" in self.symbols) #!!!
        while text:
            found=set()
            for symbol in self.symbols:
                if replaced_symbol not in symbol:
                    continue
                if text[letind:].startswith(symbol):
                    found.add(symbol)
                    
            if found:
                
                symbs=sorted(found,key=len,reverse=True)[0]
                if symbs[-1]==replaced_symbol:
                    newtext+=new_symbol
                    text=text[lo:]    
                else:
                    ls=len(symbs)
                    newtext+=text[:ls]
                    text=text[ls:]
                    
            elif text[letind:].startswith(replaced_symbol):
                newtext+=new_symbol
                text=text[lo:]
            else:
                newtext+=text[:1]
                text=text[1:]
        return newtext

    def lexicalize(self,text,temp_ops=set(),replace_forbidden=True):
        letind=0
        #print("Checking",text)
        #input()
        end=len(text)
        while letind<end:
            found=set()
            for symbol in self.symbols:
                if text[letind:].startswith(symbol):
                    found.add(symbol)
                    
            if found:
                symb=sorted(found,key=len,reverse=True)[0]
                if symb in self.operators:
                    self.pops[symb]=self.operators[symb]
                elif symb in self.funcs:
                    self.pfuncs[symb]=self.funcs[symb]
                letind+=len(symb)
            else:
                for symbol in self.forbidden_symbols:
                    if text[letind:].startswith(symbol):
                        found.add(symbol)
                        
                        
                if found:
                    symb=sorted(found,key=len,reverse=True)[0]
                    newsymb=self.make_unique(symb)
                    end+=2*len(self.modifier)
                    text=text[:letind]+newsymb+text[letind+len(symb):]
                    if newsymb in self.operators:
                        self.pops[newsymb]=self.operators[newsymb]
                    elif newsymb in self.funcs:
                        self.pfuncs[newsymb]=self.funcs[newsymb]
                    letind+=len(newsymb)
                else:
                    letter=text[letind]
                    if letter not in {"."," ",",","|"} and not letter.isdigit():
                        if self.implicit_conversion:
                                if letter.isalpha() or not self.restrict_symb:
                                    self.add_polvar(letter,temp_ops=temp_ops)
                                    self.pops[letter]=self.operators[letter]
                                    if self.verbose:
                                        print(f"Converting {letter} to variable")
                                else:  
                                    raise ValueError(f"Can't implicitly convert the non-alphabetic symbol {letter} to variable. "
                                "To override this, set restrict_symb to 'False'")
                        else:
                            raise ValueError("Found illegal letter ",letter,"which is not in",self.symbols)
                            #return True,letter
                    letind+=1
        if "," in text:
            try:
                self.pfuncs[self.def_func]=self.funcs[self.def_func]
                self.pfuncs[self.def_func2]=self.funcs[self.def_func2]
            except KeyError:
                pass
        return text

    def fill_implicit(self,text):
        lenop=len(self.implicit_op)
        for lbrack,rbrack in self.brackets:
            while True:
                ot=text
                text=re.sub(f"(\d)\s+(\d)",lambda s: f"{s[1]}{self.implicit_op}{s[2]}",text)
                if ot==text:
                    break
            
            text=re.sub(f"(\d)\s*{re.escape(lbrack)}",lambda s: s[1]+f"{self.implicit_op}{lbrack}",text)
            text=re.sub(f"{re.escape(rbrack)}\s*(\d)",lambda s: f"{rbrack}{self.implicit_op}"+s[1],text)
            
            for symb,op in self.pops.items():
                if op.arity==0: 
                    #text=re.sub(f"({re.escape(symb)})\s*{re.escape(lbrack)}",lambda s: s[1]+f"{self.implicit_op}{lbrack}",text)
                    text=re.sub(f"{re.escape(rbrack)}\s*({re.escape(symb)})",lambda s: f"{rbrack}{self.implicit_op}"+s[1],text) #komplet nahrazeni, zde asi nehrozi nebezpeci
                    text=re.sub(f"({re.escape(symb)})\s*{re.escape(lbrack)}",lambda s: s[1]+f"{lbrack}",text) # Jen odstrani mezeru
                    lpad=0 
                    for ind in self.find_all(text,f"{symb}{lbrack}"):                        
                        if not self.is_substring(text,ind+lpad,len(f"{symb}"),let_pass=self.opsymbols)[0]:
                            ns=f"{symb}{self.implicit_op}{lbrack}"
                            text=text[:ind+lpad]+ns+text[ind+lpad+len(ns)-1:]
                            #text=re.sub(f"({symb2}){name}",lambda s: f"{s[1]}{self.implicit_op}{name}",text)
                            lpad+=lenop
                for symb,op in self.pfuncs.items():
                    #text=re.sub(f"({re.escape(symb)})\s*{re.escape(lbrack)}",lambda s: s[1]+f"{self.implicit_op}{lbrack}",text)
                    text=re.sub(f"{re.escape(rbrack)}\s*({re.escape(symb)})",lambda s: f"{rbrack}{self.implicit_op}"+s[1],text)
        
        for lbrack in (brack[0] for brack in self.brackets):
            for rbrack in (brack[1] for brack in self.brackets):
                text=re.sub(f"{re.escape(rbrack)}\s*{re.escape(lbrack)}",f"{rbrack}{self.implicit_op}{lbrack}",text)
        for name in self.pfuncs:
            text=re.sub(f"(\d)\s*{name}",lambda s: f"{s[1]}{self.implicit_op}{name}",text)
            for symb2,op2 in self.pops.items():
                if op2.arity==0: 
                    lpad=0 
                    for ind in self.find_all(text,f"{symb2}{name}"):      
                        if not self.is_substring(text,ind+lpad,len(f"{symb2}{name}"))[0]:
                            ns=f"{symb2}{self.implicit_op}{name}"
                            text=text[:ind+lpad]+ns+text[ind+lpad+len(ns)-1:]
                            lpad+=lenop
                    nt=None
                    while True:
                        ot=text
                        text=re.sub(f"({symb2})\s+{name}",lambda s: f"{s[1]}{self.implicit_op}{name}",text)
                        if ot==text:
                            break
        for symb,op in self.pops.items():
            if op.arity==0 or (op.arity==1 and op.acts_on=="left"):
                if op.arity==0:
                    while True:
                        #print(text,symb)
                        res=re.findall(f"({symb}(\d+))",text)
                        #res=re.findall(f"({symb}(\d+))",text.encode("unicode_escape"))
                        if res:
                           # print(res)
                            for pres in res:
                            #print(res,res[0],res[1])
                                if pres[0] not in self.pops:
                                    text=re.sub(pres[0],f"{symb}{self.implicit_op}{pres[1]}",text)
                        break
                    text=re.sub(f"(\d)\s*{symb}",lambda s: f"{s[1]}{self.implicit_op}{symb}",text)
                    #text=re.sub(f"{symb}\s*(\d)",lambda s: f"{symb}{self.implicit_op}{s[1]}",text)
                for symb2,op2 in self.pops.items():
                    if op2.arity==0 or (op2.arity==1 and op2.acts_on=="right"):
                        lpad=0
                        for ind in self.find_all(text,f"{symb}{symb2}",hungry=True):
                                if not self.is_substring(text,ind+lpad,len(f"{symb}{symb2}"))[0]:
                                    ns=f"{symb}{self.implicit_op}{symb2}"
                                    text=text[:ind+lpad]+ns+text[ind+lpad+len(ns)-1:]
                                    lpad+=lenop
                        while True:
                            ot=text
                            text=re.sub(f"{symb}\s+({symb2})",lambda s: f"{symb}{self.implicit_op}{s[1]}",text)
                            if ot==text:
                                break
            elif op.arity==1:
                if op.acts_on=="right":
                    text=re.sub(f"(\d)\s*{symb}",lambda s: f"{s[1]}{self.implicit_op}{symb}",text)
                if op.acts_on=="left":
                    text=re.sub(f"{symb}\s*(\d)",lambda s: f"{symb}{self.implicit_op}{s[1]}",text)
        for symb in self.pfuncs:
            text=re.sub(f"{symb}\{self.implicit_op}",symb,text)
        return text

    def check_save(self,text):
        ind=text.find("<--")
        if ind!=-1:
            varname=text[:ind].strip()
            rest=text[ind+3:]
            return varname,rest
        else:
            return None,text

    def substitute(self,text,eval_dict={},muleval_dict={},overwritten_funcs=set(),overwritten_ops=set(),temp_ops=set(),**kwargs):
        text,found_brackets=self.strip_brackets(text,aggresive=False) 
        multisub={}
        vars=[el for (lb,rb) in self.brackets for el in 
              re.findall(f"([Vv]ar{re.escape(lb)}(.*?){re.escape(rb)})",text)
        ]
        if vars:
            for expr,var in vars:
                for lb,rb in self.brackets:
                    if lb in var or rb in var:
                        raise ValueError("The argument of var can't contain any parentheses!")
                name=self.add_polvar(var,overwritten_funcs,overwritten_ops,temp_ops,pretty=True)
                #print(text,expr,name)
                #input()
                text=text.replace(expr,name)
        masked_text=self.mask_text_im(text)
        #var_ind=masked_text.find("|")
        subsymbs=("|","where")
        var_inds=[]
        for symb in subsymbs:
            vind=masked_text.find(symb)
            if vind!=-1:
                var_inds.append((vind,len(symb)))
        #var_ind=masked_text.find("|")
        #if var_ind!=-1:
        #var_enum=text[var_ind+1:].split(";")
        #text=text[:var_ind
        if var_inds:
            var_ind=min(var_inds)
            #var_enum=text[var_ind[0]+var_ind[1]:].split(";")
           # print("ahoj")
            var_enum=re.split(r";|\b and \b",text[var_ind[0]+var_ind[1]:].strip(";"))
            #print(var_enum)
            #print(var_enum)
            text=text[:var_ind[0]]
            #sub_list=[]
            
            for sub in (s for s in var_enum if s):
                #print("sub",sub)
                ind=sub.find(":=") 
                #assert ind!=-1, f"Wrong substitution in {sub}"
                if ind!=-1:
                    #print("Found :=",sub)
            
                    symbol,value=sub[:ind].strip(),sub[ind+2:].strip()
                    #print(symbol,value)
                 
                    
                    if "(" in symbol:
                       # print("in symbol")
                    #if value.startswith("\\"):
                        
                            #value=eval("lambda "+value[2:])
                        fname,rest=symbol.split("(")
                        fname=fname.strip()
                        
                        if not fname:
                            raise ValueError("There must be a symbol before the parentheses in function definition!")
                            
                        fargs=tuple(farg.strip() for farg in rest[:-1].split(","))
                        if any (farg.strip()=="" for farg in fargs):
                            raise ValueError("There must be no  empty argument in the function definition!")
                        nargs=len(fargs)
                        #print("||".join(map(str,(fname,rest,fargs,value))))
                        #print(value+"|"+";".join(f"{farg}:={arg}" for arg,farg in zip((4,),fargs)))
                        
                        #except SyntaxError as e:
                        #raise ValueError("lambda "+new_func[1:])
                        #    raise ValueError(str(e)+" Wrongly specified operator!")
                        
                        if fname in self.funcsymbols:
                            overwritten_funcs.add(self.funcs[fname])
                        elif fname in self.opsymbols:
                            overwritten_ops.add(self.operators[fname])
                        #print("printing",value)
                        
                        #if self.is_number(value):
                        #    value=self.convert_to_num(value)
                        #else:
                        #    tempkwargs={**kwargs}
                        #    tempkwargs["draw"]=False
                        #    value=self.parse(value,**tempkwargs)
                        #print("printing",value)
                        
                        #if isinstance(value,ImplicitVector): #tady se dá ověřit, jestli to je implicit vector. Pokud ano, tak se dá výraz postupně dosadit a vyhodnotit pro každou složku (?)
                        mt=self.mask_text_im(value)
                        
                        if ","  in mt:
                            #print("ahooj",self.brackets[0][0],value)
                            newtexts=[]
                            #print("brackeets")
                            for ind,pos in enumerate(res.start() for res in re.finditer(",",mt+",")):
                                
                                comp=value[:pos]
                                value=value[pos+1:]
                                print(pos,comp,value,"ahoj")
                            #for ind,comp in enumerate(value):
                            #    print("comp:",comp,str(comp))
                                new_fname=f"#{fname}{ind}#"
                                temp_ops.add(new_fname)
                                new_func=lambda *args,expr=comp,fargs=fargs: self.parse(expr+"|"+";".join(f"{farg}:={arg}" for arg,farg in zip(args,fargs)))
                                #new_func=lambda *args,value=value,fargs=fargs: value.eval(**{farg:arg for arg,farg in zip(args,fargs)})
                                self.add_function(new_fname,nargs,new_func,overwrite=True)

                                #newtexts.append(self.sparse(text,verbose=False))
                                newtexts.append(self.replace_unknown(fname,new_fname,text))
                            
                            text=",".join(newtexts)
   
                        else:
                            #print("nobrack")
                            temp_ops.add(fname)
                            new_func=lambda *args,value=value,fargs=fargs: self.parse(value+"|"+";".join(f"{farg}:={arg}" for arg,farg in zip(args,fargs)))
                            #new_func=lambda *args,value=value,fargs=fargs: value.eval(**{farg:arg for arg,farg in zip(args,fargs)})
                            self.add_function(fname,nargs,new_func,overwrite=True)
                        
                    else:
                        temp_ops.add(symbol)
                        if self.is_number(value):
                            value=self.convert_to_num(value)
                        else:
                            tempkwargs={**kwargs}
                            tempkwargs["draw"]=False
                            value=self.parse(value,**tempkwargs)
                            
                        if symbol in self.funcsymbols:  
                            overwritten_funcs.add(self.funcs[symbol])
                        elif symbol in self.opsymbols:
                            overwritten_ops.add(self.operators[symbol])
                        
                        if isinstance(value,ImplicitVector): #tady se dá ověřit, jestli to je implicit vector. Pokud ano, tak se dá výraz postupně dosadit a vyhodnotit pro každou složku (?)
                            newtexts=[]
                            #print("OText:",text)

                            for ind,comp in enumerate(value):
                                new_symb=f"#{symbol}{ind}#"
                                #print("new_symb",symbol,new_symb)
                                #self.add_operator(symbol,10**4,0,lambda x=comp: x,overwrite=True)    
                                self.add_operator(new_symb,10**4,0,lambda x=comp: x,overwrite=True)    
                                #print(newtexts)
                                #newtexts.append(self.sparse(text,verbose=False))
                                newtexts.append(self.replace_unknown(symbol,new_symb,text))

                            text=",".join(newtexts)
                            #print("NText:",text)
                        else:
                            self.add_operator(symbol,10**4,0,lambda value=value: value,overwrite=True)
                        
                        self.add_operator(symbol,10**4,0,lambda value=value: value,overwrite=True)
               

                else:
                    ind=sub.find("->") 
                    if ind!=-1:
                        #print("Found =",sub)
                        varname,exprs=sub[:ind].strip(),sub[ind+2:].strip()
                        brack_inds=self._filter_brackets(exprs)
                        commas=list(self.find_all(self.mask_text(exprs,brack_inds),","))
                        if commas:
                           # print("COMMAS!",list(commas))                                
                            expr_list=[]
                            last_pos=0
                            for comma_pos in commas:
                                expr_list.append(exprs[last_pos:comma_pos])
                                last_pos=comma_pos+1
                            expr_list.append(exprs[last_pos:])
                            multisub[varname]=expr_list
                            #new_texts=(re.sub(f"{varname}",f"{expr}" ,text) for expr in expr_list) #odkomentovat pro původní nested verzi
                            #text="("+",".join(new_texts)+")"
                        else:
                            text=re.sub(f"{varname}",f"{exprs}" ,text)  
                
                    
                    else:
                        sub=sub.strip()
                        #if sub.startswith("{") and sub.endswith("}"):
                        #    sub=sub.lstrip("{").rstrip("}")
                        #    mlist=sub.split(",")
                            
                        ind=sub.find("=") 

                        if ind!=-1:
                            

                            #print("Found ->",sub)
                            symbol=sub[:ind].strip()
                            if symbol.isalpha() or not self.restrict_symb:

                                self.add_polvar(symbol,overwritten_funcs,overwritten_ops,temp_ops)
                                tempkwargs={**kwargs}
                                tempkwargs["draw"]=False

                                exprs=sub[ind+1:].strip()
                                """ # Původní verze než jsme to dělali přes implicit vector 
                                brack_inds=self._filter_brackets(exprs)
                                commas=list(self.find_all(self.mask_text(exprs,brack_inds),","))
                                if False:#commas:
                                    #print("COMMAS!",list(commas))                                
                                    expr_list=[]
                                    last_pos=0
                                    for comma_pos in commas:
                                        expr_list.append(exprs[last_pos:comma_pos])
                                        last_pos=comma_pos+1
                                    expr_list.append(exprs[last_pos:])
                                    muleval_dict[symbol]=[self.parse(expr,**tempkwargs) for expr in expr_list]
                        
                                else:
                                """
                                res=self.parse(exprs,**tempkwargs)
                                if isinstance(res,ImplicitVector):
                                    muleval_dict[symbol]=[comp for comp in res]
                                else:
                                    eval_dict[symbol]=self.parse(exprs,**tempkwargs)
                                
                            else:
                                raise ValueError(f"Can't use {symbol} as a variable name, use letters of the alphabet."
                                "To override this, set restrict_symb to 'False'")
                            #print(self.eval_dict)
                        else:        
                            symbol=sub.strip()
                            
                            self.add_polvar(symbol,overwritten_funcs,overwritten_ops,temp_ops)
            if multisub: #if False: # pro nested verzi
        
                vals=([(key,comp) for comp in val] for key,val in multisub.items())
                combs=it.product(*vals)
                textcol=[]
                for subs in combs:
                    nt=text
                    for varname,expr in subs:
                        nt=re.sub(f"{varname}",f"{expr}" ,nt) #zakomentovat
                    textcol.append(nt)
                text=",".join(txt for txt in textcol)
        #print(found_brackets)
        if found_brackets:
            return "("+text+")"
        else:
            return text 
    def add_monvar(self,symbol,temp_ops=set()):
        temp_ops.add(symbol)
        self.add_operator(symbol,10**7,0,
        lambda value=Monom(1,{symbol:1}):value)
    def add_polvar(self,symbol,overwritten_funcs=set(),overwritten_ops=set(),temp_ops=set(),pretty=False):
        
        if pretty:
            var=PrettyVariable(symbol)
        else:
            var=Variable(symbol)
        symbol=str(var) 
        temp_ops.add(symbol)
        if symbol in self.funcsymbols: 
            overwritten_funcs.add(self.funcs[symbol])
        elif symbol in self.symbols:
            overwritten_ops.add(self.operators[symbol])
        self.add_operator(symbol,10**7,0,lambda value=var:value,overwrite=True)    
        return symbol

    def check_brackets(self,text):
        
        #bracknums=[[0,0] for _ in self.brackets]
        last_left=[]
        for letind,letter in enumerate(text):
            for ind,brack in enumerate(self.brackets):
                if letter==brack[0]:
                  #  bracknums[ind][0]+=1
                    last_left.append(ind)
                elif  letter==brack[1]:
                   # bracknums[ind][1]+=1
                    if not last_left:# bracknums[ind][1]>bracknums[ind][0]:
                        raise ValueError(f"Spatne uzavrene zavorky (pravych zavorek je do pozice {letind+1} vic nez levych)!")
                    if ind!=last_left.pop():
                        raise ValueError(f"Zkrizene zavorky na pozici {letind+1}!")
                
                
        if  last_left:
        #if any(num[1]!=num[0] for num in bracknums):
            raise ValueError("Spatne uzavrene zavorky (levych je vic nez pravych)!")
    


    def preprocess(self,text):
        #text=re.sub("(--)+","+",text) #sudy pocet minusu -> jedno plus
        for op in self.pops:#operators:
            if self.operators[op].arity<2:
                continue
            for op2 in self.pops:#operators:
                if self.operators[op2].arity<2:
                    continue
                symb=self.unmake_unique_text(op)
                symb2=self.unmake_unique_text(op2)
                #print(op,op2)
                
                if symb in text and symb2 in text:
                    if self.operators[op].pri>self.operators[op2].pri:
                        text=re.sub(f"{re.escape(symb)}{re.escape(symb2)}(-?\d+\.?\d*)",lambda s: f"{symb}({symb2}{s[1]})",text) #*- -> *(-)
        
        return text 

    def parse(self,text:str,draw=False,verbose=None,implicit_conversion=None,preparse=None,
              implicit_ops=None,force_func_brackets=None,fractions=None,solve_w_complex=None,
              solformat=None,replace_forbidden=True,**kwargs):

        #___________________setup_____________________________
        eval_dict={}
        meval_dict={}
        opops,opfuncs,self.pops,self.pfuncs=self.pops,self.pfuncs,{},{}
        #for_dict={}
        overwritten_funcs=set()
        overwritten_ops=set()
        self.partial_results=[]
        self.orig_text=text
        temp_ops=set()
        orig_impl=self.implicit_conversion
        if implicit_conversion is not None:
            self.implicit_conversion=implicit_conversion
        orig_verbose=self.verbose
        if verbose is not None:
            self.verbose=verbose
        orig_preparse=self.preparse
        if preparse is not None:
            self.preparse=preparse
        orig_implicit_ops=self.implicit_ops
        if implicit_ops is not None:
            self.implicit_ops=implicit_ops
        orig_force_func_brackets=self.force_func_brackets
        if force_func_brackets is not None:
            self.force_func_brackets=force_func_brackets
        orig_fractions=self.fractions
        if fractions is not None:
            self.fractions=fractions
        orig_sol_w_complex=self.solve_w_complex
        if solve_w_complex is not None:
            self.solve_w_complex=solve_w_complex
        orig_solformat=self.solformat
        if solformat is not None:
            self.solformat=solformat
        if fractions is not None:
            self.fractions=fractions
        #__________________budeme ukladat vyraz do symbolu?______________
        savevar,text=self.check_save(text) 
        #__________________substituce | _________________
        text=self.substitute(text,eval_dict=eval_dict,muleval_dict=meval_dict,
            	            temp_ops=temp_ops,overwritten_funcs=overwritten_funcs,overwritten_ops=overwritten_ops,replace_forbidden=replace_forbidden,**kwargs)
        #__________________preparsing a kontrola zavorek_________________ 
        if self.preparse:
            
            text=re.sub("(--)+","+",text) #sudy pocet minusu -> jedno plus
            
            self.check_brackets(text)

        #__________________nahrada zakazanych symbolu, aby nekolidovali operatory a numerika (treba -7 je cislo a zaroven operator - aplikovany na 7) -> "-" se zmeni na "#-#"_________________
        if False:#replace_forbidden:
            
            for symb in self.forbidden_symbols:
                if symb in text:
                    newsymb=self.make_unique(symb)
                    if not newsymb in text:
                        if newsymb in self.symbols:    
                            text=text.replace(symb,newsymb)

        
        #__________________nahrazeni neznamych pismenek za operatory ci funkce_________________
        text=self.lexicalize(text,temp_ops,replace_forbidden=replace_forbidden)
        #__________________dalsi preparsing, typu *-   -> *(-)
        if self.preparse:
            
            text=self.preprocess(text)
            
        try:    
                            #if found_illegal:
                                    #    raise ValueError(f"The expression {text} contains an undefined symbol '{illegal_symbol}'!")

            #__________________doplni implicitni symbol, vetsinou krat, takze z 5x se stane 5*x
            if self.implicit_ops:
                self.pops[self.implicit_op]=self.operators[self.implicit_op]
                text=self.fill_implicit(text)
            
            self.pops=dict(sorted(self.pops.items(),key=self._sort_op_key))
            self.pfuncs=dict(sorted(self.pfuncs.items(),key=lambda item:-len(item[0])))
            
            if self.preparse:
                if self.verbose:
                    
                    print("Preparsed: ", self.unmake_unique_text(text))
            res=self._parse(text)
            
                
            if eval_dict:
                if isinstance(res,GeneralObject):
                    res=res.eval(**eval_dict)
                    
            if meval_dict:
                try:
                    vals=( [(key,comp) for comp in val] for key,val in meval_dict.items())
                    combs=list(it.product(*vals))
                    args=[dict(comb) for comb in combs]
                    if isinstance(res,GeneralObject):    
                        res=ImplicitVector(*(res.eval(**arg) for arg in args))
                    else:
                        res=ImplicitVector(*(res for arg in args))
                except TypeError as e:
                    print(str(e))
                    pass
        except OverflowError:
            raise ValueError("The result is too large")       
        finally:
            #print("finally")
            for op in set(temp_ops):
                #print("op:",op)
                try:
                    if op in self.opsymbols:
                        #print("op op",op)
                        self.opsymbols.remove(op)
                        del self.operators[op]
                    elif op in self.funcsymbols:
                        #print("func:",op)
                        self.funcsymbols.remove(op)
                        del self.funcs[op]
                    self.symbols.remove(op)
                    
                except (KeyError,IndexError):
                    print("error!",op)
                    pass 
            for func in overwritten_funcs:
                    self.funcs[func.symbol]=func
                    self.funcsymbols.add(func.symbol)
                    self.add_symbol(func.symbol)
            for op in overwritten_ops:
                    self.operators[op.symbol]=op
                    self.opsymbols.add(op.symbol)
                    self.add_symbol(op.symbol)
            self._sort_ops()
            self.implicit_conversion=orig_impl
            self.verbose=orig_verbose
            self.preparse=orig_preparse
            self.implicit_ops=orig_implicit_ops
            self.fractions=orig_fractions
            self.solformat=orig_solformat
            self.sol_w_complex=orig_sol_w_complex
            self.force_func_brackets=orig_force_func_brackets
            self.pops,self.pfuncs=opops,opfuncs
            self.__cur_error_depth=0
        #print(eval_dict)
        if savevar is not None:
            #print("adding",savevar)
            if "(" in savevar:
                fname,rest=savevar.split("(")
                fname=fname.strip()
                if not fname:
                    raise ValueError("There must be a symbol before the parentheses in function definition!")
                            
                fargs=tuple(farg.strip() for farg in rest[:-1].split(","))
                if any (farg.strip()=="" for farg in fargs):
                    raise ValueError("There must be no  empty argument in the function definition!")
                nargs=len(fargs)
                #new_func=lambda *args,value=res,fargs=fargs: self.parse(value+"|"+";".join(f"{farg}:={arg}" for arg,farg in zip(args,fargs)))
                new_func=lambda *args,value=res,fargs=fargs: value.eval(**{farg:arg for arg,farg in zip(args,fargs)})
                self.add_function(fname,nargs,new_func,overwrite=True)
            else:
                self.add_operator(savevar,10**7,0,lambda value=res:value,overwrite=True)    
        if draw:
            self.draw_parse_tree(**kwargs)
        if False:#isinstance(res,GeneralObject):
            if res.contains_const():
                cvars=res.get_const_vars()
                res=res.eval(**{c:0 for c in cvars})
        return res


    def nicepower(self,expr,lsymb="<sup>",rsymb="</sup>"):
        expr=re.sub("\^(-?\d+)",lambda s: f"<sup>{s[1]}</sup>",expr)
        ind=expr.find("^")      
        while ind!=-1:
            
            for bind,(lp,rp) in enumerate(self.brackets):
                #if any(b in lp or b in rp for b in ("{","}")): #tohle může vést k nečekanému chování!
                    
                if expr[ind+1]==lp:
                    loc=self._filter_brackets(expr[ind+1:])[bind][0][1]
                    expr=expr[:ind]+"<sup>"+expr[ind+1:][1:loc]+"</sup>"+expr[ind+2:][loc:]
                    break
            else:
                break
            ind=expr.find("^")
        expr=re.sub("<sup>",lsymb,expr)
        expr=re.sub("</sup>",rsymb,expr)
        return expr

    def sparse(self,*args,html=False,**kwargs):
        res=self.parse(*args,**kwargs)
        #res=res.eval()
        if isinstance(res,(int,float)):
            return gen_obj.truncate_number(res)
        else:
            expr=str(res)
            if html:
                return self.nicepower(expr)
            else:
                return expr
        
    def _parse(self,text:str,level=(0,0),lbuf=None,rbuf=None):
        #print(f"_parser got {text}")
        self.partial_results.append([level,self.unmake_unique_text(text),"sym",None])
        #print(text)
        #print(text)
        text,found_brackets=self.strip_brackets(text,aggresive=True)
        #print(text)
        if self.is_final(text): 
        #    print(f"final {text}")
            #return(str(float(text)))
            #print(lbuf,rbuf)
            if lbuf is not None or rbuf is not None:
                
                if text[0] not in ("+","-"):
                    #print(text,lbuf,rbuf)

                    raise ValueError("Spatne zadany vyraz! (Nejspis unarni operator mezi dvema argumenty nebo unarni operator na spatne strane)")
            res=self.convert_to_final(text) # a zde by se mohli vracet obecnejsi algebraicke objekty
            #print(text,res)
            self.partial_results.append([level,res,"res",None])
            return res


        #self.partial_results.append((level,self.unmake_unique_text(text),"sym"))
        
        brack_inds=self._filter_brackets(text)
        masked_text=self.mask_text(text,brack_inds)
        
        if masked_text.find("|")!=-1 or masked_text.find("where")!=-1 :
            #print("masked_text:",masked_text,".......",text)  
            #return self.parse(text,replace_forbidden=False)
            return self.parse(text,replace_forbidden=True)
        elif masked_text.find(",")!=-1 and masked_text.find("..")==-1:
            if found_brackets:
                text=f"{self.def_func2}({text})"
            else:
                text=f"{self.def_func}({text})"
        
        else:
            for symbol,op in self.pops.items():
               
                for ind in self.find_all(masked_text,symbol,direction=op.direction):
                  
                    is_conflict,_=self.is_substring(text,ind,len(symbol),let_pass={symbol})
                    if is_conflict: 
                        continue

                    #print("Passed conflict check: ",symbol)
                    lp=text[:ind]
                    rp=text[ind+len(symbol):]
                    #print(symbol)
                    #print(f"Before eval: lp = {lp}, rp= {rp}")
                    
                    if op.arity==0:
                        res=op.intr()
                        #print(f"res {res} , {type(res)}")
                        self.partial_results.append([level,res,"res",None])
                        return res
                    else:
                        
                        if op.arity==1:
                            res=self.expr_eval1(lp,rp,op,level=level,lbuf=lbuf,rbuf=rbuf)
                            self.partial_results.append([level,res,"res",None])
                            return res
                    #       print(lp,rp)
                        
                        elif op.arity==2:
                            #print("ar2",op.symbol,lp,rp,lbuf,rbuf)
                            if lp=="" and lbuf==None and op.unary_overload:
                            #    print("Overloading")
                            #       print("Overloading right")
                                    if op.unary_overload.acts_on=="left":                            
                                        raise ValueError(f"Wrong direction of action! The operator {symbol} acts to the left.")
                                    else:
                                        res=self.expr_eval1(lp,rp,op.unary_overload,level=level,lbuf=lbuf,rbuf=rbuf)
                                        self.partial_results.append([level,res,"res",None])
                                        return res
                            elif rp=="" and rbuf==None and op.unary_overload:
                                    #print(lbuf,rbuf,lp,rp)
                                    if op.unary_overload.acts_on=="right":
                                        raise ValueError(f"Wrong direction of action! The operator {symbol} acts to the right.")
                                    else:
                                        res=self.expr_eval1(lp,rp,op.unary_overload,level=level,lbuf=lbuf,rbuf=rbuf)
                                        self.partial_results.append([level,res,"res",None])
                                        return res
                            else:
                                res=self.expr_eval2(lp,rp,op,level=level,lbuf=lbuf,rbuf=rbuf)
                                self.partial_results.append([level,res,"res",None])
                                return res
        #print(text)
        res=self.eval_funcs(text,level=level)
        self.partial_results.append([level,res,"res",None])
        return res

    def expr_eval1(self,lp,rp,op,level=(0,0),lbuf=None,rbuf=None):
            #print(f"Unary, expr= {expr},op= {op.symbol}")
            #if self.is_number(expr):
            #    return op.intr(float(expr))
            #else:
            #return self.expr_eval1(self._parse(expr),op)
        xl,yl=level
        if op.acts_on=="left": 
            if not rp:
                if not lp:
                    if lbuf is not None:
                        return op.intr(lbuf)
                    else:
                        print(lp,rp,op.symbol,lbuf,rbuf,sep="|")
                        raise ValueError("Wrong input!")
                return op.intr(self._parse(lp,level=(xl,yl+1),lbuf=lbuf,rbuf=rbuf))
                                    
            else:
                #print(lp,rp,lbuf,rbuf)
                res=op.intr(self._parse(lp,level=(xl,yl+1),lbuf=lbuf,rbuf=None))
                #print("res",res)
                #assert lbuf==!None,f"No lbuf in {lp}!"
                    
                #print(lp,rp)
                return self._parse(rp,level=(xl+1/2**yl,yl+1),lbuf=res,rbuf=rbuf)
        else: #if op.acts_on=="right"
            if not lp:
                if not rp:
                    if rbuf is not None:
                        return op.intr(rbuf)
                    else:
                        print(lp,rp,op.symbol,lbuf,rbuf,sep="|")
                        raise ValueError("Wrong input!")
                return op.intr(self._parse(rp,level=(xl,yl+1),lbuf=lbuf,rbuf=rbuf))
            else: #if lp
                #assert lbuf==None,f"Conflicting lbuf in {lp}!"
                #print(rp,None,rbuf)
                res=op.intr(self._parse(rp,level=(xl,yl+1),lbuf=None,rbuf=rbuf))
                #print(lp,rp)
                return self._parse(lp,level=(xl-1/2**yl,yl+1),lbuf=lbuf,rbuf=res)

    def expr_eval2(self,lp,rp,op,level=(0,0),lbuf=None,rbuf=None):
        xl,yl=level
        
        if lp=="":
            
            if lbuf is not None:
                self.partial_results[-1][3]=self.unmake_unique_text(op.symbol)    
                return op.intr(lbuf,self._parse(rp,level=(xl+1/2**yl,yl+1),rbuf=rbuf))   
            elif self.fill_id and op.identity:
                lp=str(op.id_rep)
            else:
                raise ValueError("Wrong input!",lp,rp)
        elif rp=="":
            if rbuf is not None:
                self.partial_results[-1][3]=self.unmake_unique_text(op.symbol+str(rbuf))    
                return op.intr(self._parse(lp,level=(xl-1/2**yl,yl+1),lbuf=lbuf),rbuf)
            elif self.fill_id and op.identity:
                rp=str(op.id_rep)
            else:
                raise ValueError("Wrong input!")
        #return self.expr_eval2(self._parse(lp),self._parse(rp),op)
        self.partial_results[-1][3]=self.unmake_unique_text(op.symbol)    
        return op.intr(self._parse(lp,level=(xl-1/2**yl,yl+1),lbuf=lbuf),self._parse(rp,level=(xl+1/2**yl,yl+1),rbuf=rbuf))

    def eval_funcs(self,text,level=(0,0)):
        #print(f"eval_funcs got {text}")
        #for symb in self.funcs:
        for symb in self.pfuncs:
            if text.startswith(symb):
                
                if self.force_func_brackets:
                    if not any((text[len(symb):].strip()).startswith(lbra) for lbra,_ in self.brackets):
                        raise ValueError("The functional arguments are not enclosed by bracekts! Set force_func_brackets to False if you wish to ignore this")
                #print(f"Stripping brackets from: {text}")
                text,_=self.strip_brackets(text[len(symb):],aggresive=False)
                #if self.funcs[symb]._textual:
                #    return self.funcs[symb].intr(text)
                #print(text[len(symb):])
                #print(f"Stripped brackets: {text}")

                brack_inds=self._filter_brackets(text)
                #print(text,brack_inds)
                commas=self.find_all(self.mask_text(text,brack_inds),",")
                #print(f"Masking text and getting comma indices: {self.mask_text(text,brack_inds)}")
                expr_list=[]
                last_pos=0
                for comma_pos in commas:
                    expr_list.append(text[last_pos:comma_pos])
                    last_pos=comma_pos+1

                expr_list.append(text[last_pos:])
                arity=self.funcs[symb].arity
                arg_num=len(expr_list)
                
                if arity != "*":
                    if arg_num!=arity:
                        raise ValueError("Spatny pocet argumentu!") 
                xl,yl=level
                if arg_num%2==0:
                    xlevels=list(range(-(arg_num//2),0))+list(range(1,arg_num//2+1))
                else:
                    xlevels=list(range(-(arg_num//2),arg_num//2+1))
                return self.funcs[symb].intr(*(self._parse(el,level=(xl+xlevels[ind]/arg_num**yl,yl+1)) 
                                            for ind,el in enumerate(expr_list)))
        raise ValueError("Wrong input (are the functional arguments bracekted correctly?)")
        #print(f"eval_funcs is returning {text}")
        

    def is_final(self,text):
        if any(text.startswith(pl) for pl in self.pollist):
            return True
        try:
            float(text)
        except ValueError:
            
            return False
        else:
            return True
    
    def convert_to_final(self,text):
        for pl in self.pollist:
            if text==pl:
                return SPolynomial(pl,0,1)
        try:
            res=int(text)
        except ValueError:
            res=float(text)
        return res

    def simple_tree(self,part_res):
        part_res=sorted(part_res,key=lambda t: t[0][1])
        tree=[]#list((el[2],el[1]part_res)
        levels=[]
        
        for el in part_res:
            level=el[0][1]
            if level not in levels:
                levels.append(level)
                tree.append([])
            if el[2]=="sym":
                tree[level].append({"sym":el})
                for rel in part_res:
                    if rel[0]==el[0] and rel[2]=="res":
                        tree[level][-1]["res"]=rel
                        break
        #print(expr_by_levels)
        return tree
        
        

    def print_simple_tree(self,indtensity=5,row_height=1):
        indtensity=indtensity
        row_height=row_height
        base=3
        tree=self.simple_tree(self.partial_results)
        for exprs in tree:
            #print(exprs)
            #input()
            sorted_expr=sorted(exprs,key=lambda el: el["res"][0][1] )
            #print(sorted_expr)
            cp=0
            
            for el in sorted_expr:
                x=el["sym"][0][0]
                y=el["sym"][0][1]
                #pos=(base+x)
                pos=(base+round(x))
                true_space=(pos)*indtensity
            #    print(" "*(true_space-cp),end="")
                text="{0} [={1},:{2}]".format(el["sym"][1],el["res"][1],el["sym"][3] if el["sym"][3] else "")
             #   print(text,end="")
                cp=true_space+len(text)
                
                
            print("\n"*row_height,end="")
        
    def draw_parse_tree(self,pr=None,xc=0.5,yc=0.9,size0=20,adapt_size=False,zepr=None,save=False,filename=None):
        if pr is None:
            tree=self.simple_tree(self.partial_results)
        else:
            tree=self.simple_tree(pr)
        if save:
            matplotlib.use("Agg")
            #x0,y0=0.15,1.1
        else:
            matplotlib.use(self.backend)
            #matplotlib.use(self.backend)
            #x0,y0=0.2,1.1
        #print(matplotlib.get_backend())
        height=0.2
        width=0.2
        size_exp=1.1
        fig,ax=plt.subplots(figsize=(15,7))
        mng = plt.get_current_fig_manager()
        try:
            mng.window.state("zoomed")
        except AttributeError:
            pass
        ax.text(0.5,1.1,self.orig_text,size=size0,color="green",horizontalalignment="center",verticalalignment="center")

        for exprs in tree:
            
            sorted_expr=sorted(exprs,key=lambda el: el["res"][0][1] )
            for el in sorted_expr:
                level=el["sym"][0][1]
                x=xc+el["sym"][0][0]*width
                y=yc-level*height
                text1=el["sym"][1]
                res=el["res"][1]
                if type(el["res"][1])==float:
                    text2=str(res) if type(res)==int else gen_obj.truncate_number(res,3)
                else:
                    text2=str(res) 
                text2="$"+self.nicepower(text2,lsymb="^{",rsymb="}")+"$"
                
                if adapt_size:
                    size=size0*size_exp**(-level)
                else:
                    size=size0
                ax.text(x,y,text1,size=size,horizontalalignment="center",verticalalignment="center")
                ax.text(x,y+height/2,text2,size=size*0.8,color="red",horizontalalignment="center",verticalalignment="center")
                if el["sym"][3]:
                    x=xc+(el["sym"][0][0])*width
                    y=yc-(level+1)*height
                    ax.text(x,y,el["sym"][3],size=size,color="blue",horizontalalignment="center",verticalalignment="center")
            for spine in ax.spines.values():
                spine.set_visible(False)
            ax.xaxis.set_ticks([])        
            ax.yaxis.set_ticks([]) 
            level+=1
        if save:
            plt.savefig(filename+str(self.plot_ind)+".png")
            plt.close()
            self.plot_ind+=1
            matplotlib.use(self.backend)
        else:
            if self.backend=="Qt5Agg":
                plt.ion()
            plt.show()
                    

        #plt.close(fig)

    def draw_exprs(self,*expressions,verbose=False,**kwargs):
        for ind,expr in enumerate(expressions):
            self.parse(expr,draw=True,save=True,filename=f"SavedResults/Expr. # {ind+1}",**kwargs,verbose=verbose)

    def map_iter(self,func,init_point,n=100):
        try:
            if n<0 or int(n)!=n:
                raise ValueError("The third argument must be a non-negative natural number!")
        except TypeError:
            raise ValueError("The third argument must be a non-negative natural number!")

        if not isinstance(func,GeneralObject):
            func=Monom(func)
        if not isinstance(func,Vector):
            try:
                func=Vector(*func)
            except TypeError:
                func=Vector(func)
        if not isinstance(init_point,Vector):
            try:
                init_point=Vector(*init_point)
            except TypeError:
                init_point=Vector(init_point)
        if func.dim!=init_point.dim:
            raise ValueError("The dimensions of the mapping function and of the initial point must be equal!")
        p=init_point
        res=[p]
        vars=sorted(func.get_vars())
        for _ in range(n):
            p=func.eval(**{var:comp for var,comp in zip(vars,p)})
            res.append(p)
        return Vector(*res)

    def multi_map_iter(self,funcs,init_points,n=100):
        res_array=[]
        try:
            iter(funcs)
        except TypeError:
            funcs=[funcs]
        try:
            iter(init_points)
        except TypeError:
            init_points=[init_points]

        for func in funcs:
            for init_point in init_points:
                res_array.append(self.map_iter(func,init_point,n))
        return res_array
    
    """
    def dif_eq(self,func,init_point,n=100):
        from scipy.integrate import odeint
        try:
            if n<0 or int(n)!=n:
                raise ValueError("The third argument must be a non-negative natural number!")
        except TypeError:
            raise ValueError("The third argument must be a non-negative natural number!")

        if not isinstance(func,GeneralObject):
            func=Monom(func)
        if not isinstance(func,Vector):
            try:
                func=Vector(*func)
            except TypeError:
                func=Vector(func)
        if not isinstance(init_point,Vector):
            try:
                init_point=Vector(*init_point)
            except TypeError:
                init_point=Vector(init_point)
        if func.dim!=init_point.dim:
            raise ValueError("The dimensions of the mapping function and of the initial point must be equal!")
        tarr=np.linspace(0,1.0,n)
        adjvars=func.get_vars().difference({"t"})
        ordvars={var for var in adjvars if not var.startswith("d")}
        dvars=adjvars.difference(ordvars)
        for var in dvars:
            if var[1:] not in ordvars:
                ordvars.add(var[1:])

        if not adjvars:
            adjvars={"y"}
        vars=sorted(ordvars)+sorted(dvars)
        
        #print(vars)
        def f(y,t):
            return [*func.eval(**{var:comp for var,comp in zip(["t"]+vars,[t]+list(y))})]
        return np.c_[tarr,odeint(f,init_point,tarr,args=tuple())]
        #return Vector(*res)

    def multi_dif_eq(self,funcs,init_points,n=100):
        res_array=[]
        try:
            iter(funcs)
        except TypeError:
            funcs=[funcs]
        try:
            iter(init_points)
        except TypeError:
            init_points=[init_points]

        for func in funcs:
            for init_point in init_points:
                print(func)
                print(init_point)
                res_array.append(self.dif_eq(func,init_point,n))
        return res_array

    def pairwise_dif_eq(self,args,n=100):
        #args=list(args)
        #print(args)
        try:
            funcs,init_points=zip(*args)
        except TypeError:
            funcs,init_points=args
            funcs=[funcs]
            init_points=[init_points]
        
        names=[f"f= {func},init: {init}" for func,init in zip(funcs,init_points)]
        res_array=[]
        for func,init_point in zip(funcs,init_points):
            res_array.append(self.dif_eq(func,init_point,n))
        return names,res_array
    """

    def load_standard(self):
        def power(x,y):
            res=x**y
            if isinstance(res,complex):
                return Complex(res)
            else:
                return res
        def eq_op(lhs,rhs):
            return Equation(lhs,rhs,self.solve)
        def eqsyst_op(left,right):
            return SystemOfEqs(left,right,system_solver=self.solve_linear_system)
        def solve_func(eq):
            if not isinstance(eq,(Equation,SystemOfEqs)):
                raise ValueError("The argument must be an equation or a system of equations")
            else:
                return eq.eval(fractions=self.fractions)

        #def ineq_op(lhs,rhs,symbol="<"):
        #    return self.solve_ineq(lhs,rhs,symbol=symbol)
        def ineq_op(lhs,rhs,symbol="<"):
            return Inequality(lhs,rhs,self.solve_ineq,symbol=symbol)
        self.def_func="IV"
        self.def_func2="V"
        self.add_operator("*",100,2,lambda x,y: x*y,identity={1})
        self.add_operator("/",110,2,lambda x,y: x/y,direction="left")
        #self.add_operator("=",-10,2,self.solve)
        self.add_operator("=",-10,2,eq_op)
        for symb in ("<",">","<=",">="):
            self.add_operator(symb,-10,2,ft.partial(ineq_op,symbol=symb))
        self.add_operator(";",-11,2,eqsyst_op)
        self.add_function("solve",1,solve_func)
        self.add_operator("==",-10,2,lambda x,y: "True" if x==y else "False")
        self.add_operator("+",10,2,lambda x,y: x+y,identity={0},unary_func=(lambda x: x,"right"))
        self.add_operator("-",12,2,lambda x,y: x-y,identity={0},direction="left",unary_func=(lambda x: -x,"right"))
        self.add_operator("#",-1,1,lambda x: x**2,acts_on="right")
        self.add_operator("^",200,2,power)
        self.add_operator("**",200,2,power)
        self.add_operator("e",1000,0,lambda: math.e)
        self.add_operator("pi",1000,0,lambda: math.pi) 
        self.add_operator("dx",1000,0,lambda: Variable("dx")) 
        self.add_operator("dy",1000,0,lambda: Variable("dy")) 
        self.add_operator("π",1000,0,lambda: math.pi)
        self.add_operator("i",1000,0,lambda: Complex(0,1))
        
        for card in ("♦","♣","♥","♠"):
            self.add_operator(card,1000,0,lambda s=card: Variable(s))    

        
        def ldiv(x,y): 
            if type(x) in [int,float,Monom,Quotient,Variable]:
                x=Polynomial(x)
            try:
                return x.longdiv(y)[0]
            except AttributeError:
                raise AttributeError("Both arguments must be polynomials!")
        def ldiv_modulo(x,y):
            if type(x) in [int,float,Monom,Quotient,Variable]:
                x=Polynomial(x)
            elif isinstance(x,Vector):
                return x.__class__(*(ldiv_modulo(comp,y) for comp in x))
            try:
                return x.longdiv(y)[1]
            except AttributeError:
                raise AttributeError("Both arguments must be polynomials!")
        def ldivcomb(x,y):
            if type(x) in [int,float,Monom,Quotient,Variable]:
                x=Polynomial(x)
            try:
                return Vector(*x.longdiv(y)[:2])
            except AttributeError:
                raise AttributeError("Both arguments must be polynomials!")                 

        def smart_deriv(obj,def_var=None,full_set=None):
            if isinstance(obj,(int,float)):
                return 0
            vars=obj.get_vars()
            
            if len(vars)==1:
                
                return obj.deriv(next(iter(vars)))
            else:
                return ImplicitVector(*obj.grad())
            
        def grad(obj):
            if isinstance(obj,(int,float)):
                return Vector(0)
            else:
                return obj.grad()
        

        def deriv(obj,*variables):
            for variable in variables:
                if isinstance(obj,(int,float)):
                    return 0
                if not isinstance(variable,Variable):
                    variable=Variable(repr(variable))
                obj=obj.deriv(variable)
            return obj
        def simplify(obj):
            if not isinstance(obj,GeneralObject):
                return obj
            else:
                return obj.simplify()
        def dot(v,w):
            
            try:
                return gen_obj.Dot_product(v,w).eval()
                # v.dot(w)
            except AttributeError:
                raise ValueError("The left argument of the function 'dot' is not a vector!") #Tohle se asi teď už může smazat, ale radši to nechme, kdyby něco.
            
        def cross(v,w):
            try:
                return gen_obj.Cross_product(v,w).eval()
            except AttributeError:
                raise ValueError("The left argument of the function 'cross' is not a vector!")
            
        def curl(v):
            try:
                #print(v,type(v))
                return gen_obj.Curl(v).eval()
            except AttributeError:
                raise ValueError("The argument of the function 'curl' is not a vector!")
        def veclap(v):
            return gen_obj.VectorLaplace(v).eval()

        def lap(obj):
            if isinstance(obj,(int,float)):
                return 0
            else:
                return obj.laplace()
        
        def divergence(v):
            return gen_obj.Divergence(v).eval()
            #if not isinstance(v,Vector):
            #    raise ValueError("The argument of the function 'div' must be a vector!")
            #return v.divergence()
        def matmul(obj1,obj2):
            try:
                #print(obj1,obj2)
                return gen_obj.Matmul(obj1,obj2).eval()
                #return dot(obj1,obj2.T())
            except AttributeError:
                raise TypeError("The arguments must be vectors!")
        def FT(v):
                return gen_obj.FourierTransform(v).eval()
        def IFT(v):
                return gen_obj.InverseFourierTransform(v).eval()
        def norm(obj):
            if isinstance(obj,(Complex,Vector)):
                return obj.norm()
            elif isinstance(obj,(int,float)):
                return abs(obj)
            else:
                return gen_obj.Abs(obj).eval()
        
        def sum_obj(obj):
            #try:
            return gen_obj.Sum(obj).eval()
            #except TypeError:
            #    return obj
        def mul_obj(obj):
            return gen_obj.Mul(obj).eval()
            #try:
            #    return ft.reduce(mul,obj)
            #except TypeError:
            #    return obj
        def trans(obj):
            try:
                return gen_obj.Transpose(obj).eval()
            except AttributeError:
                
                return obj
        def fact_decomp(*nums):
            def _fact_decomp(num):
                if isinstance(num,Vector):
                    return num.__class__(*(_fact_decomp(comp) for comp in num))
                if isinstance(num,(Variable,Monom,Polynomial)):
                    vars=num.get_vars()
                    if len(vars)<1:
                        num=float(num)
                        num=int(num) if int(num)==num else num
                    elif len(vars)>1:
                        raise NotImplementedError("Can't decompose multivariate polynomials yet")
                    else:
                        use_complex=False
                    #  print(vars)
                        var=next(iter(vars))
                        x=Variable(var)
                        roots=self.solve(num,0,format=False,fractions=False)
                        
                        if not isinstance(roots,Vector):
                            roots=[roots]
                        if use_complex:
                            return Vector(*(x-root for root in roots))
                        else:
                            real_roots=[root.real for root in roots if (not isinstance(root,Complex)) or root.isreal()]
                            compl_roots=[root for root in roots if root not in real_roots]
                            
                            partners=[]
                            found=set()
                            for root in compl_roots:
                                if root not in found:
                                    for root2 in compl_roots:
                                        if root is not root2 and root2 not in found:
                                            if root==root2.conj() or root==-root2.conj():
                                                
                                                partners.append((root,root2))
                                                found.add(root)
                                                found.add(root2)
                            #print(*(repr(el) for el in real_roots))
                            #print("rp",real_roots,compl_roots,partners)
                            return Vector(*(x-root for root in real_roots),*( ((x-r1)*(x-r2)).real.simplify() for (r1,r2) in partners))
                if not isinstance(num,int):
                    raise TypeError("The argument must be a positive integer or a polynomial in one variable!")
                return Vector(*prime_factor.prime_fact(num))

            if len(nums)==1:
                return _fact_decomp(nums[0])
            else:
                return ImplicitVector(*(_fact_decomp(num) for num in nums))   

        def mulparent_candidates(*nums):
            def _mulparent_candidates(num,lens="all"):
                if isinstance(num,Vector):
                    return num.__class__(*(Vector(*_mulparent_candidates(comp,lens=lens)) for comp in num))
                if isinstance(lens,Vector):
                    return lens.__class__(*(Vector(*_mulparent_candidates(num,lens=comp)) for comp in lens))
                factors=(group for group in fact_decomp(num))
                return ImplicitVector(*(Vector(*el) for el in prime_factor.sortedmulgroups(list(factors)) if lens=="all" or len(el)==lens))

            if len(nums)==1:
                return _mulparent_candidates(nums[0])
            elif len(nums)==2:
                return _mulparent_candidates(nums[0],lens=nums[1])
                #return ImplicitVector(*(Vector(*_mulparent_candidates(num)) for num in nums))
            else:
                raise ValueError("The function 'parents' may have at most two arguments!")

        def choose(n1,n2):
            res=gen_obj.Choose(n1,n2).eval()
            if isinstance(res,float):
                if int(res)==res:
                    res=int(res)
            return res
        def factorial(num):
            #try:
            #     return math.factorial(num)        
            #except ValueError,TypeError:
            return gen_obj.Gamma(num+1).eval()
        def gamma(num):
            #try:
            #     return math.factorial(num)        
            #except ValueError,TypeError:
            return gen_obj.Gamma(num).eval()
        def stirling(n):
            return (gen_obj.Sqrt(2*np.pi*n)*(n/np.e)**n).eval()
        def mag(x):
            return 2/3*(gen_obj.Log(x)-9.1).eval()
        def M0(y):
            return 10**(3/2*y+9.1)
        def poisson(k,λ):
            return gen_obj.Poisson(k,λ).eval()
        def beta(num1,num2):
            return gen_obj.Beta(num1,num2).eval()#(gen_obj.Gamma(num1)*gen_obj.Gamma(num2)/(gen_obj.Gamma(num1+num2))).eval()
        def normal(x,μ=0,σ=1):
            return gen_obj.NormalDist(x,μ,σ).eval()
        def lognormal(x,μ=0,σ=1):
            return gen_obj.LogNormal(x,μ,σ).eval()
        def binom(p,n,k):
            return gen_obj.Binomial(p,n,k).eval()
        def maxwell_boltzmann(v,T=273.15,m=28.97):
            #π=3.14159265359
            #σ=gen_obj.Sqrt(1000*8.31446261815324*T/m).eval()
            return gen_obj.Maxwell_Boltzmann(v,T,m).eval()# 1/(2*π*σ**2)*v**2*normal(v,0,σ) 
        def random_sample(n=1,x0=0,x1=-1,func=1):
            if n==0:
                n=1
            return gen_obj.RandomSample(n,x0,x1,func).eval()# 1/(2*π*σ**2)*v**2*normal(v,0,σ) 
        def random_normal_sample(n=1,μ=0,σ=1):
            if n==0:
                n=1
            return gen_obj.RandomSample(n,0,-1,gen_obj.NormalDist(Variable("s"),μ,σ)).eval()# 1/(2*π*σ**2)*v**2*normal(v,0,σ)
        def phase(z):
            z=Complex(z)
            return gen_obj.Angle(z.real,z.imag).eval()
        def flat(v):
            try:
                return v.flat()
            except AttributeError:
                return v
        def take_real(x):
            try:
                if isinstance(x,GeneralObject):
                    #print("isinst")
                    return x.simplify().real
                return x.real
            except AttributeError:
                return x
        def take_imag(x):
            try:
                if isinstance(x,GeneralObject):
                    return x.simplify().imag
                return x.imag
            except AttributeError:
                return Monom(0)
        def get_component(v,ind):
            return v[ind]
        def common_factor(x):
            x=Polynomial(x)
            return Vector(*x.factor_variables())
        def gcd(*x):
            if len(x)<2:
                raise ValueError("The number of arguments must be greater than 1!")
            x1,*x=x
            while len(x)>0:
                x2,*x=x
                x1=gen_obj.Gcd(x1,x2).eval()
            return x1
        def lcm(*x):
            if len(x)<2:
                raise ValueError("The number of arguments must be greater than 1!")
            x1,*x=x
            while len(x)>0:
                x2,*x=x
                x1=gen_obj.Lcm(x1,x2).eval()
            return x1
        def id(n):
            return Vector(*(Vector(*(1 if i==k else 0 for i in range(n))) for k in range(n) ))
        def fmax(*pars):
            return (ft.reduce(gen_obj.Max,pars)).eval()
        def fmin(*pars):
            return (ft.reduce(gen_obj.Min,pars)).eval()
        
                
        self.add_operator(":",-5,2,ldivcomb)
        self.add_operator("%",-4,2,ldiv_modulo)
        self.add_operator("'",193,1,smart_deriv,acts_on="left")
        def extended(func,ext):
            def ext_func(x):
                if isinstance(x,(int,float)): 
                    return func(x)
                elif isinstance(x,GeneralObject):
                    return ext(x).eval()
            return ext_func
        ln=extended(math.log,gen_obj.Ln)

        
        self.add_operator("!",194,1,factorial,acts_on="left")
        self.add_function("sin",1,extended(math.sin,gen_obj.Sin))
        self.add_function("cos",1,extended(math.cos,gen_obj.Cos))
        self.add_function("asin",1,extended(math.asin,gen_obj.Asin))
        self.add_function("acos",1,extended(math.acos,gen_obj.Acos))
        self.add_function("sqrt",1,lambda x: gen_obj.Sqrt(x).eval())
        self.add_function("√",1,lambda x: gen_obj.Sqrt(x).eval())
        self.add_function("exp",1,extended(math.exp,gen_obj.Exp))
        self.add_function("ln",1,ln)
        self.add_function("log",1,extended(math.log10,gen_obj.Log))
        self.add_function("log_",2,lambda base,x: ln(x)/ln(base))
        self.add_function("tan",1,extended(math.tan,gen_obj.Tan))
        self.add_function("atan",1,extended(math.atan,gen_obj.Atan))
        self.add_function("deg",1,lambda x: x*180/math.pi)
        self.add_function("rad",1,lambda x: x*math.pi/180)
        self.add_function("abs",1,extended(abs,gen_obj.Abs))
        self.add_function("angle",2,lambda x,y: gen_obj.Angle(x,y).eval())
        self.add_function("phase",1,phase)
        self.add_function("conj",1,lambda z: gen_obj.Conj(z).eval())
        self.add_function("J1",2,lambda n,z: gen_obj.Bessel(z,n).eval())
        self.add_function("J2",2,lambda n,z: gen_obj.SecondBessel(z,n).eval())
        self.add_function("E1",2,lambda phi,z: gen_obj.Ellip1(phi,z).eval())
        self.add_function("E2",2,lambda phi,z: gen_obj.Ellip2(phi,z).eval())
        self.add_function("CE1",1,lambda z: gen_obj.CompEllip1(z).eval())
        self.add_function("CE2",1,lambda z: gen_obj.CompEllip2(z).eval())
        self.add_function("sign",1,lambda x: gen_obj.Sign(x).eval())
        self.add_function("one",1,lambda x: gen_obj.Const(x,1))
        self.add_function("null",1,lambda x: gen_obj.Const(x,0))
        self.add_function("const",2,lambda x,y: gen_obj.Const(x,y))
        self.add_function("normal","*",normal)
        self.add_function("lognormal","*",lognormal)
        self.add_function("binom",3,binom)
        self.add_function("gammadist",3,lambda x,α,β: gen_obj.GammaDist(x,α,β).eval())
        self.add_function("betadist",3,lambda x,α,β: gen_obj.BetaDist(x,α,β).eval())
        self.add_function("gammafunc",1,gamma)
        self.add_function("stirling",1,stirling)
        self.add_function("mag",1,mag)
        self.add_function("moment",1,M0)
        self.add_function("betafunc",2,beta)
        self.add_function("poisson",2,poisson)
        self.add_function("MB","*",maxwell_boltzmann)
        self.add_function("rand","*",random_sample)
        self.add_function("nrand","*",random_normal_sample)
        self.add_function("BT",3,lambda x,y,z: gen_obj.BayesianFormula(x,y,z).eval())
        self.add_function("H",1,lambda x: gen_obj.Heaviside(x).eval())
        self.add_function("V","*",Vector)
        self.add_function("IV","*",ImplicitVector)
        self.add_function("dot",2,dot)
        self.add_function("cross",2,cross)
        self.add_function("Q","*",Quotient)
        self.add_function("grad",1,grad)
        self.add_function("FT",1,FT) 
        self.add_function("IFT",1,IFT) 
        self.add_operator("×",101,2,cross)
        self.add_operator("·",101,2,dot)
        self.add_operator("∇",105,1,grad,acts_on="right")
        self.add_operator("∇·",101,1,divergence,acts_on="right")
        self.add_operator("∇×",102,1,curl,acts_on="right")
        self.add_operator("Δ",104,1,lap,acts_on="right")
        self.add_operator("_",200,2,get_component,direction="left")
        self.add_operator("..",-20,2,lambda x,y: gen_obj.VectorGenerator(x,y).eval(),direction="right")
        self.add_operator("...",-20,2,lambda x,y: gen_obj.VectorGenerator(x,y).eval(),direction="right")
        self.add_function("Id",1,id)
        self.add_function("list",3,lambda x,y,z: gen_obj.VectorGenerator(gen_obj.Vector(x,y),z).eval())
        self.add_function("sinh",1,extended(math.sinh,gen_obj.Sinh))
        self.add_function("cosh",1,extended(math.cosh,gen_obj.Cosh))
        self.add_function("tanh",1,extended(math.tanh,gen_obj.Tanh))
        self.add_function("asinh",1,extended(math.asinh,gen_obj.Asinh))
        self.add_function("acosh",1,extended(math.acosh,gen_obj.Acosh))
        self.add_function("atanh",1,extended(math.atanh,gen_obj.Atanh))
        self.add_function("erf",1,lambda x: gen_obj.Erf(x).eval())
        self.add_function("erfc",1,lambda x: gen_obj.Erfc(x).eval())
        self.add_function("curl",1,curl)
        self.add_function("rot",1,curl)
        self.add_function("div",1,divergence)
        self.add_function("lap",1,lap)
        self.add_function("veclap",1,veclap)
        self.add_function("matmul",2,matmul)
        self.add_function("der","*",deriv)
        self.add_function("simp",1,simplify)
        self.add_function("trans",1,trans)
        self.add_function("decomp","*",fact_decomp)
        self.add_function("factor",1,common_factor)
        self.add_function("parents","*",mulparent_candidates)
        self.add_function("gcd","*",gcd)
        self.add_function("lcm","*",lcm)
        self.add_function("norm",1,norm)
        self.add_function("sum",1,sum_obj)
        self.add_function("mul",1,mul_obj)
        self.add_function("flat",1,flat)
        self.add_function("choose",2,choose)
        self.add_function("max","*",fmax)
        self.add_function("min","*",fmin)

        self.add_function("real",1,take_real)
        self.add_function("Real",1,lambda x: gen_obj.Real(x).eval())
        self.add_function("imag",1,take_imag)
        self.add_function("ode","*",de.DifEq)
        self.add_function("bvp","*",de.BVP)
        self.add_function("ivp","*",de.DifEq)
    
        self.add_function("map","*",lambda *args: de.DifEq(*args,is_map=True))
        self.add_function("bifurc","*",lambda *args: de.DifEq(*args,is_map=True,bifurc=True))
        #def map_iter_plot(points,)

        def lin_interface(*args):
            lhs_list,rhs_list=[],[]
            if args:
                try:
                    for i in range(0,len(args),2):
                        lhs_list.append(args[i])
                        rhs_list.append(args[i+1])
                except IndexError:
                    raise ValueError("Must supply an even number of arguments!")
            else:
                return 0
            return self.solve_linear_system(lhs_list,rhs_list)
        self.add_function("solve_syst","*",lin_interface)
                    
        self.add_brackets("()")
        self.add_brackets("{}")
        self.add_brackets("[]")

    def solve(self,lhs,rhs,format=None,fractions=None):
        
        if format is None:
            format=self.solformat
        if fractions is None:
            fractions=self.fractions
        if isinstance(lhs,Quotient):
            lhs,rhs=lhs.p,rhs*lhs.q
        if isinstance(rhs,Quotient):
            lhs,rhs=lhs*rhs.q,rhs.p

        lhs,rhs=lhs-rhs,0
        if isinstance(lhs,Complex):
            if lhs.isreal():
                lhs=lhs.real
            else:
                raise NotImplementedError("Can't solve equations with complex coefficients yet")
        lhs=Polynomial(lhs)
        #print("lhs...",lhs,repr(lhs))
        if lhs.order>4:
            return appr_polsolve.absolve(lhs)
        else:            
            #print(len(lvar),len(rvar))
            #if len(lvar)==len(rvar)==0:
            if (not lhs.varset):
                
                if lhs!=0:
                    return "Never true"
                else:
                    return "Always true"
            else:
                #print("1")
                #print(lhs,lhs.varset)
                lvar=next(iter(lhs.varset))
                if len(lhs.varset)>1:
                    if lhs.order>-1:
                        raise NotImplementedError("Rovnice musi byt ve stejnych promennych, nelinearni soustavy zatim neumime resit.")
                    else:
                        return self.solve_linear_system([lhs],[0],format=format,fractions=fractions)

            #print(lvar)
            a0,a1,a2,a3,a4=lhs.get_coefs(lvar,at_least=4)#.values()
            #print("coefs:",a0,a1,a2,a3,a4)
           # print(a1,b1,c1)
           # print(a2,b2,c2)
            if a4==0:
                if a3==0:
                    if a2==0:
                        return solver.solve_linear(a1,a0,fractions=fractions)
                    else:
                        return solver.solve_quadratic(a2,a1,a0,format=format,fractions=fractions,use_complex=self.solve_w_complex)
                else:
                    return solver.solve_cubic(a3,a2,a1,a0,fractions=fractions,use_complex=self.solve_w_complex)
            else:
                return solver.solve_quartic(a4,a3,a2,a1,a0,fractions=fractions,use_complex=self.solve_w_complex)

    def solve_ineq(self,lhs,rhs,symbol="<",format=None,fractions=None):
            roots=Equation(lhs,rhs,self.solve).eval(format=format,fractions=False)
            dif=lhs-rhs
            sign=1
            #fon=fractions or (fractions is not None and self.fractions)
            dif=Polynomial(dif)
            vars=dif.get_vars()
            if len(vars)>1:
                raise NotImplementedError("Inequalities in multiple variables haven't been implemented yet.")
            if symbol in (">",">="):
                sign=-1
            if roots =="Always true" or roots == "Never True":
                
                vars=dif.get_vars()
                var="x" if not vars else next(iter(vars))
                if roots=="Always true":
                    if symbol in ["<=",">="]:
                        return f"{var} ∈ (-∞,∞)"
                    else:
                        return f"{var} ∈ {{}}"
                elif roots=="Never true":
                    return f"{var} ∈ (-∞,∞)" if dif.eval(**{var:0})*sign<0 else f"{var} ∈ {{}}"
            try:
                s=0
                print(roots)
                #if fon:
                #    rroots=sorted([root.real for root in roots if (not isinstance(root,Complex)) or root.isreal()],key=lambda x: x.real.eval())
                #else:
                rroots=sorted({root.real for root in roots if (not isinstance(root,Complex)) or root.isreal()})
            except TypeError as e:
                print(s,e)
                rroots=[roots]
            bounds=["-∞"]+rroots+["∞"]
            print(bounds)
            intervals=[(bounds[i],bounds[i+1]) for i,_ in enumerate(bounds[:-1]) ]
            lb,rb="(",")"
            if symbol in ("<=",">="):
                lb,rb="[","]"
            tr=[]
            for interval in intervals:
                l,r=interval
                lpar,rpar=lb,rb
                if l=="-∞":
                    lpar="("
                    if r=="∞":
                        rpar=")"
                        val=0
                    else:
                        val=r-1
                elif r=="∞":
                    rpar=")"
                    val=l+1
                else:
                    val=0.5*(l+r)
                var=next(iter(dif.get_vars()))
                if isinstance(val,Quotient):
                    val=val.eval()
                
                if dif.eval(**{var:val})*sign<0:   
                    
                    if l not in ("-∞","∞"):
                        #if fon:
                        #    l=str(l)
                        #else:
                        l=gen_obj.truncate_number(l)
                    if r not in ("-∞","∞"):
                        #if fon:
                        #r=str(r)
                        #else:
                        r=gen_obj.truncate_number(r)
                    
                    tr.append(lpar+l+","+r+rpar)
            #return tr
            if not tr:
                tr="{}",
            return f"{var} ∈ "+"∪".join(tr)
    def extract_bounds(self,result):
        result=result[result.find("∈")+1:].strip()
        if result=="{}":
            return [("empty","empty")],False
        elif result=="(-∞,∞)":
            return [("-∞","∞")],False
        else:
            closed= ("[" in result or "]" in result)
            ints=result.split("∪")
            bounds=[]
            for int in ints:
                b1,b2=int.split(",")
                b1=b1[1:]
                b2=b2[:-1]
                if b1!="-∞":
                    b1=float(b1)
                if b2!="∞":
                    b2=float(b2)
                bounds.append((b1,b2))
            return bounds,closed
    def solve_linear_system(self,lhs_list,rhs_list,format=None,fractions=None):
        proc_lhs_list=[]
        tol=1.e-5
        Arows=[]
        b=[]
        varset=set()
        for ind,(lhs,rhs) in enumerate(zip(lhs_list,rhs_list)):
            fractions=self.fractions
            if isinstance(lhs,Quotient):
                lhs,rhs=lhs.p,rhs*lhs.q
            if isinstance(rhs,Quotient):
                lhs,rhs=lhs*rhs.q,rhs.p
            lhs,rhs=lhs-rhs,0
            if isinstance(lhs,Complex):
                if lhs.isreal():
                    lhs=lhs.real
                else:
                    raise NotImplementedError("Can't solve equations with complex coefficients yet")
            lhs=Polynomial(lhs)
            if lhs.order>1:
                raise NotImplementedError("Can't solve non-linerar systems yet")
            proc_lhs_list.append(lhs)
            varset.update(lhs.varset)
        for lhs in proc_lhs_list:
       #     print("varset",varset)
            coefs,const=lhs.get_linear_coefs(sorted(list(varset)))
            Arows.append(list(coefs.values()))
            b.append(-const)
        A=np.array(Arows);b=np.array(b)
        #print(A,b)
        solution,solvable,kernel=QRLU.QRsolve(A,b)
        solution=solution[:,0]
        solvable=solvable[0]
        kernel=np.transpose(kernel)
        for i,k in enumerate(kernel):
            for j,_ in enumerate(k):
                if k[j]!=0:
                    kernel[i]=kernel[i]/k[j]
                    
                    break
        for i,sol in enumerate(solution):
            if abs(sol-int(sol))<tol:
                    solution[i]=int(sol)
        for i,k in enumerate(kernel):
            for j,kj in enumerate(k):
                if abs(kj-int(kj))<tol:
                        kernel[i][j]=int(kj)
        return Solution(solution,kernel,solvable,varset)
        """
        return Vector(Vector(Variable("Solution:") if solvable  else Variable("Not exactly solvable,"
        " showing approximate (LS) solution:"),Vector(*solution)),
        Vector(Variable("Kernel:"),Vector(*(Vector(*k) for k in np.transpose(kernel)))))
        """ 


    def longdiv(self,s1,s2=None,draw=False,implicit_conversion=True,**kwargs):
        if s2==None:
            #print(s1,"|",s2)

            lind=s1.find("|")
            if lind==-1:
                lind=len(s1)
            col_num=len([el for el in self.find_all(s1,":") if el<lind])
            
            if col_num!=1:
                raise ValueError("The expression must contain exactly  one ':' sign!")
            else:
                dind=s1.find(":")
                s1,s2=s1[:dind],s1[dind+1:]
             #   print(s1,"|",s2)
                if any(t.strip()=="" for t in (s1,s2)):
                    raise ValueError("There have to non-empty expressions on each side of the division sign!")
                
                if "|" in s2:
                    _,sub,*_=s2.split("|")
                    s1+=" | "+sub

        pol1=self.parse(s1,implicit_conversion=implicit_conversion,verbose=self.verbose,**kwargs)
        pol2=self.parse(s2,implicit_conversion=implicit_conversion,verbose=self.verbose,**kwargs)
        if isinstance(pol1,(int,float,Variable,Monom,Quotient)):
            pol1=Polynomial(pol1)
        if isinstance(pol2,(int,float,Variable,Monom,Quotient)):
            pol2=Polynomial(pol2)
        
        res,zbytek,partial=pol1.longdiv(pol2)
        #print(partial)
        spartial=list(map(lambda t: tuple(map(str,t)),partial))
        if draw:
            self.draw_ldiv(spartial,str(pol2),**kwargs)
        return str(res),str(zbytek),list(map(lambda t: tuple(map(str,t)),partial))
    
    def draw_ldiv(self,partial,pol2,use_colors=True,color_mode=1,save=False,filename=None,size=18):
        assert len(partial)>0,"Chybi mezivysledky!"
    
        if save:
            matplotlib.use("Agg")
            x0,y0=0.15,1.1
        else:
            matplotlib.use(self.backend)
            #matplotlib.use(self.backend)
            x0,y0=0.2,1.1
  
        fig,ax=plt.subplots(figsize=(15,7))
        mng = plt.get_current_fig_manager()
        try:
            mng.window.state("zoomed")
        except AttributeError:
            pass
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.xaxis.set_ticks([])        
        ax.yaxis.set_ticks([]) 
        def subpow(text): return re.sub("\^\(?(-?\d+)\)?",lambda s: f"${{}}^{{{s[1]}}}$",text)
        nice_partial=[(subpow(text1),subpow(text2)) for text1,text2 in partial ]
        #print(nice_partial)
        #ax.text(0.5,1.1,self.orig_text,size=size0,color="green",horizontalalignment="center",verticalalignment="center")
        font = {'family' : 'monospace',
        'weight' : 'normal',
        'size'   : size}
        matplotlib.rc('font', **font)
        #print(partial)
        nice_lhs=nice_partial[0][0]
        adjl=len(list(self.find_all(nice_lhs,"^")))
        adjconst=7
        if use_colors:
            colors=["blue","green","orange","darkmagenta"]*100
        else:
            colors=["black"]*400   
        #xlevel0=len("{0} : {1} = ".format(lhs,pol2))
        pol2=subpow(str(pol2))
        adjm=len(list(self.find_all(pol2,"^")))
        xlevel0=len(" : {0} = ".format(pol2))-adjm*7
        xlevel=xlevel0
        respos=0
        toplevel=0
        ylevel=toplevel
        height=0.1
        #size=15
        #print("{0} : {1} = {2}".format(lhs,pol2,rhs))
        if color_mode==1:
            lcolor=colors[0]
        elif color_mode==2:
            lcolor="black"
        ax.text(x0,y0+ylevel*height,"{0}".format(nice_lhs),size=size,color=lcolor,
                horizontalalignment="right",verticalalignment="center")
        ax.text(x0,y0+ylevel*height," : {0} = ".format(pol2),size=size,color="black",
                horizontalalignment="left",verticalalignment="center")
        ylevel+=1
        
        last_color=colors[0]
        if color_mode==1:
            colors=colors[1:]
        elif color_mode==2:
            last_color=colors[0]
            colors=colors
        if len(partial)==1:
            ax.text(x0,y0-toplevel*height," "*xlevel+str(0),size=size,color="black",
                    horizontalalignment="left",verticalalignment="center") 
        for (r,s),color in zip(nice_partial[1:],colors):
            if color_mode==1:
                lcolor=color
                rcolor=last_color
            elif color_mode==2:
                lcolor=color
                rcolor=color
            #print(r,s,"{0:>{ljust}}".format(r,ljust=ljust)) "x{0:>{ljust}}".format(r,ljust=ljust)
            adjl=len(list(self.find_all(r,"^")))
            ax.text(x0,y0-ylevel*height,"{0}".format(r),size=size,color=lcolor,
            horizontalalignment="right",verticalalignment="center")
            adjl=len(list(self.find_all(s,"^")))
            ax.text(x0,y0-toplevel*height," "*xlevel+" {0}".format(s[respos:]),size=size,color=rcolor,
                    horizontalalignment="left",verticalalignment="center")            
            respos=len(s)
            #print(respos)
            xlevel=xlevel0+respos-adjl*adjconst
            last_color=color
            ylevel+=1


        #if respos==0:
        #  ax.text(x0,y0+toplevel," "*xlevel+" {}".format(0),size=size,color=lcolor,horizontalalignment="center",verticalalignment="center")            
        if save:
            plt.savefig(filename+str(self.plot_ind)+".png",bbox_inches='tight')
            plt.close()
            self.plot_ind+=1
            matplotlib.use(self.backend)
                

        else:
            if self.backend=="Qt5Agg":
                plt.ion()
            plt.show()
        #plt.close(fig)
        

    def ldiv_and_save(self,*exprlist,**kwargs):
        for ind,expr in enumerate(exprlist):
            
            self.longdiv(expr,draw=True,save=True,filename=f"SavedResults/Deleni #{ind+1}",**kwargs)



    #import numpy
    def plot(self,*objs,llim=None,ulim=None,t0=0,t1=2*math.pi,tden=1000,pden=None,style="",save=False,filename=None,polar=None,
    use_grid=None,contours=None,hold=False,anim=False,pointanim=False,equal=False,realtime=False,realimag=False,v_slice=None,legend=True,scale=False,loglog=False,hlines=[],vlines=[],animt=None,**at):
        self.ac=self.AnimControl()
        show2D=False
        waveanim=False
        is_open=False
        varset=set() #for label
        wavearr=[]
        new_objs=[]
        names=[]
        wave_names=[]
        rinds=[]
        iinds=[]
        imagvals=[]
        parinds=[]
        
        if save:
            matplotlib.use("Agg")
            if filename==None:
                filename="SavedResults/Plot"
            else:
                filename="SavedResults/"+filename
        else:
           #matplotlib.use(self.backend)
           matplotlib.use(self.backend)
        try:
            x0,y0=llim
        except TypeError:
            x0=llim
            y0=None
        try:
            x1,y1=ulim
        except TypeError:
            x1=ulim
            y1=None
        newobjs=[]
        for obj in objs:
            
            if isinstance(obj,Solution):
                for sol in obj:
                    newobjs.append(sol)
            elif isinstance(obj,Inequality):
                newobjs.append(obj)
                newobjs.append(obj.lhs)
                newobjs.append(obj.rhs)
                

            elif isinstance(obj,Equation):
                newobjs.append(obj.eval(format=False,fractions=False))
            elif isinstance(obj,SystemOfEqs):
                for eq in obj:
                    newobjs.append(eq.eval(format=False,fractions=False))

            elif isinstance(obj,Vector) and obj.dim>2:
                if v_slice not in (None,"all"):
                    newobjs.append(Vector(obj[v_slice[0]],obj[v_slice[1]]))
                else:
                    if not obj.get_vars():
                        newobjs.append(obj.eval())
                    else:
                        for comb in it.combinations(obj,2):
                            newobjs.append(Vector(*comb))
            elif isinstance(obj,GeneralObject):
                newobjs.append(obj.eval())
            else:
                newobjs.append(obj)
        objs=newobjs
        
        if all(isinstance(obj,(int,float,complex))  or isinstance(obj,GeneralObject) and not obj.get_vars() for obj in objs): #konstanty (hlavně na vypisování kořenů)
            
            real,imag,names=[],[],[]
            if all(type(obj) is Vector and obj.dim==2 for obj in objs):
                draw_pairs=True
                numlist=objs
            else:
                draw_pairs=False
                numlist=Vector(*objs).flat()
            if scale:
                numlist=numlist/max(abs(numlist))
            if not any(isinstance(obj,(complex,Complex)) for obj in numlist) and (not draw_pairs):
                fig,ax=plt.subplots()
                nl=len(numlist)
                xdata=list(range(nl))
                self.plotdata.append([xdata,numlist])
                
                ax.plot(numlist,"." if not style else style,markersize=9)
                ax.set_xticks(range(0,nl+1,max(1,nl//25)))
                ax.set_xlabel("k")
                ax.set_ylabel("f(k)")
                if loglog:
                    ax.set_xscale("log")
                    ax.set_yscale("log")
                if equal:
                    plt.gca().set_aspect("equal")
                if use_grid:
                    plt.gca().grid(True)
                if save:
                    plt.savefig(filename+str(self.plot_ind)+".png")
                    self.plot_ind+=1
                    plt.close()
                    matplotlib.use(self.backend)
                else:
                #  if self.backend=="Qt5Agg":
                #      plt.ion()
                    #plt.ion()
                    
                    plt.show()
            else:
                for obj in numlist:
                    
                    #obj=Complex(obj).eval()
                    #print(obj)
                    #print(obj.real)
                    if draw_pairs:
                        real.append(obj[1]);imag.append(obj[2]);names.append(str(obj))
                    else:
                        real.append(obj.real if isinstance(obj,(complex,Complex,Vector)) else obj);imag.append(obj.imag if isinstance(obj,(complex,Complex,Vector)) else 0)
                        names.append(str(obj.regularized() if isinstance(obj,(complex,Complex,Vector)) else obj))
                origx0,origx1,origy0,origy1=x0,x1,y0,y1
                x0=min(real) if x0 is None else x0;x1=max(real) if x1 is None else x1
                
                if x0>x1:
                    x0,x1=x1,x0
                x0=x0-0.15*(x1-x0);x1=x1+0.15*(x1-x0)
                tol=Complex(None).tol
                if (x1-x0)<tol:
                    x0-=1
                    x1+=1
                if y0 is None:
                    y0=min(imag);
                if y1 is None:
                    y1=max(imag) 
                    y0=y0-0.15*abs(y1-y0);y1=y1+0.15*abs(y1-y0)
                xdif=x1-x0
                if y1-y0<1/7*xdif:
                    y0,y1=(y0+y1)/2-1/14*xdif,(y1+y0)/2+1/14*xdif
                ydif=y1-y0
                if x1-x0<1/7*ydif:
                    x0,x1=(x0+x1)/2-1/14*ydif,(x1+x0)/2+1/14*ydif
                xdif,ydif=x1-x0,y1-y0
                
                if origx0 is not None:
                    x0=origx0
                if origx1 is not None:
                    x1=origx1
                if origy0 is not None:
                    y0=origy0
                if origy1 is not None:
                    y1=origy1

                fig,ax=plt.subplots()
                def sign(x):
                    if x>0:
                        return 1
                    elif x<0:
                        return -1
                    else:
                        return 0
                
                for r,i,name in zip(real,imag,names):
                    ax.plot(r,i,"o")
                if loglog:
                    ax.set_xscale("log")
                    ax.set_yscale("log")
                if len(real)<10:
                    for r,i,name in zip(real,imag,names):
                        xtext,ytext=r-0.05*(x1-x0)*sign(r),i-0.05*(y1-y0)*sign(i)
                        if abs(xtext)<1/14*(xdif):
                            xtext+=1/14*(xdif)
                        if abs(ytext)<1/14*(ydif):
                            ytext+=1/14*(ydif)
                        print(name)
                        ax.annotate(name,(r,i),ha="center",xytext=(xtext,ytext),xycoords="data")
                
                ax.set_xlim(x0,x1);ax.set_ylim(y0,y1)
                if len(real)>=10:
                    ax.set_xticks(np.linspace(x0,x1,10));ax.set_yticks(np.linspace(y0,y1,10))
                ax.axhline(0,lw=1,color="black")
                ax.axvline(0,lw=1,color="black")
                if equal!=False:
                    plt.gca().set_aspect("equal")
                if use_grid!=False:
                    plt.gca().grid(True)
                for hline in hlines:
                    ax.axhline(hline,lw=1)
                for vline in vlines:
                    ax.axvline(vline,lw=1)
                tsfig,tsax=plt.subplots()
                
                nl=len(numlist)
                xdata=list(range(nl))
                if draw_pairs:
                    c1,c2=list(zip(*numlist))
                    self.plotdata.append([xdata,c1])
                    tsax.plot(c1,"." if not style else style,label="Comp. 1",markersize=9)
                    self.plotdata.append([xdata,c2])
                    tsax.plot(c2,"." if not style else style,label="Comp. 2",markersize=9)
                else:
                    self.plotdata.append([xdata,numlist.real])
                    tsax.plot(numlist.real,"." if not style else style,label="Real part",markersize=9)
                    self.plotdata.append([xdata,numlist.imag])
                    tsax.plot(numlist.imag,"." if not style else style,label="Imaginary part",markersize=9)
                
                tsax.set_xticks(range(0,nl+1,max(1,nl//25)))
                tsax.set_xlabel("k")
                tsax.set_ylabel("f(k)")
                if loglog:
                    tsax.set_xscale("log")
                    tsax.set_yscale("log")
                plt.legend()
                if save:
                    for fig_ind in (fig.number,tsfig.number):
                        plt.figure(fig_ind)
                        plt.savefig(filename+str(self.plot_ind)+".png")
                        self.plot_ind+=1
                        plt.close()
                    matplotlib.use(self.backend)
                else:
                #  if self.backend=="Qt5Agg":
                #      plt.ion()
                    #plt.ion()
                    tsfig.show()
                    fig.show()
                    #plt.show()
            #plt.show()
            return 

        for ind,obj in enumerate(objs): #Rozdělení potenciálně komplexních členů na reálnou a imaginární část
            if isinstance(obj,GeneralObject) and obj.contains_complex():
                
                new_objs.append(obj)
                new_objs.append(obj)
                rinds.append(ind)
                iinds.append(ind+1)
                nice_str="$"+self.nicepower(str(obj),lsymb="^{",rsymb="}")+"$"
                names.append("Real part of "+nice_str)
                names.append("Imaginary part of "+nice_str)
      
                #print(obj.real,obj.imag,objs)
            else:
                new_objs.append(obj)
                names.append("$"+self.nicepower(str(obj),lsymb="^{",rsymb="}")+"$")
        objs=new_objs
        for obj in objs:  #Automatické určení mezí pro polynomy
           
            if isinstance(obj,Polynomial) and x0==None==x1:
                try:          
                    roots=self.solve(obj,0,format=False,fractions=False)
                    try:
                        iter(roots)
                    except TypeError:
                        roots=[roots]
                    realroots=[]
                    for root in roots:
                        if isinstance(root,GeneralObject):
                            root=root.eval()
                        if isinstance(root,(int,float)): 
                            realroots.append(root)
                        if isinstance(root,Complex) and root.isreal():
                            realroots.append(root.real)
                    #print(realroots)
                    if realroots:
                        if use_grid is None:
                            #print("using grid")
                            use_grid=True
                        x0=min(realroots+[x0]) if x0 is not None else min(realroots)
                        x1=max(realroots+[x1]) if x1 is not None else max(realroots)
                        if x0==x1:
                            if x0<0:
                                if x0<-1:
                                    x1=-x0
                                else:
                                    x0,x1=-1,1
                            else:
                                if x0>1:
                                    x0,x1=-x0,x0
                                else:
                                    x0,x1=-1,1
                        else:
                            d=x1-x0
                            x0=x0-d/3
                            x1=x1+d/3
                except (ValueError,NotImplementedError):
                    pass         
        try:
            x0=x0[0]
        except TypeError:
            pass
        try:
            x1=x1[0]
        except TypeError:
            pass
        
        for ind,obj in enumerate(objs):
            
            if isinstance(obj,GeneralObject):
                if not isinstance(obj,Inequality):
                    obj=obj.eval(**at)
                
            if not isinstance(obj,GeneralObject):
                obj=Monom(obj)
            vars=obj.get_vars()
            lv=len(vars)

            if lv==0: #Nalezení jmen proměnných a vyřešení různých případů v závislosti dle jejich počtu
              #  waveanim=False
                vars="x"
                
            elif not isinstance(obj,de.DifEq) and (lv>1 and not isinstance(obj,(Equation)) or lv>2):         
                #print(type(obj))
                if not (anim and lv==2): #a co rovnice?!
                    if ind not in iinds:
                        if pden is None:       
                            ps=None
                        else:
                            ps=pden if pden<250 else 250
                        
                        self.plot2(obj,llim=llim,ulim=ulim,pden=ps,style=style,save=save,realimag=realimag,filename="2DPlot",use_grid=use_grid,polar=polar,contours=contours,equal=equal,
                        hold=True,realtime=realtime,v_slice=v_slice,anim=anim,animt=animt,t0=t0,t1=t1,**at)
                        show2D =True
                    continue
                else:
                    if lv==2:
                        waveanim=True        
                        if animt not in vars: #počítá s tím, že ve vars není nikdy None. Pokud by to byl problém, tak přidat explicitní if not animt is None
                            timevar,spacevar=sorted(vars) 
                        else:
                            timevar=animt
                            spacevar=next(iter(vars.difference({timevar,})))

                        if spacevar in ("f","fi","Φ","α","β","γ","φ","θ"):
                            if polar is None:
                                polar=True        
            else:
                #waveanim=False
                var=next(iter(vars))
                if var in ("f","fi","Φ","α","β","γ","φ","θ"):
                    if polar is None:
                        polar=True

            if not is_open:
                if polar:
                    if hold:
                        if self.backend!="TkAgg":
                            raise ValueError("The hold feature doesn't work with the current backend. Please got to ../calc_settings.json and set backend to 'TkAgg'")
                        fig=plt.gcf()
                    else:
                        fig=plt.figure()
                    ax = fig.add_subplot(111, projection='polar')
                    if x0 is None:
                        x0=0
                    if x1 is None:
                        x1=2*math.pi
                    if equal!=False:
                        plt.gca().set_aspect("equal")
                else:
                    if hold:
                        if self.backend!="TkAgg":
                            raise ValueError("This hold feature doesn't work with the current backend. Please got to ../calc_settings.json and set backend to 'TkAgg'")
                        fig,ax=plt.gcf(),plt.gca()
                    else:
                        #fig=self.fig=Figure(figsize=(5,5), dpi=100)
                        #ax=fig.add_subplot(111)
                        fig,ax=plt.subplots()
                    if all(isinstance(comp,Vector) for comp in objs):
                        if x0 is None:
                            x0=0
                        if x1 is None:
                            x1=2*math.pi

                    else:
                        if x0 is None:
                            x0=-1
                        if x1 is None:
                            x1=1
                    if equal==True:
                        plt.gca().set_aspect("equal")

                if x0>x1:
                    x0,x1=x1,x0
                is_open=True
            if pden==None:
                pden=1000
            if waveanim:
                var=spacevar #sorted(vars)[-1]
            else:
                var=sorted(vars)[0]
            varset.add(var)
            if t0 is None:
                t0=0
            if t1 is None:
                t1=2*np.pi
            if t1<t0:
                t0,t1=t1,t0
            if tden is None or tden<0:
                tden=100
            time_range=np.linspace(t0,t1,tden)
            
            grid=np.array([x0+(x1-x0)*i/pden for i in range(pden+1)]) #np.linspace(x0,x1,)
            vals=[]
            
            #print(ind,rinds,iinds)
            if isinstance(obj,de.DifEq):
                #curves,names,var_list=obj.calculate()
                curves,names,var_list=obj.calculate(t0=t0,t1=t1,n=tden)
                #print(curves,names,var_list)
                #plt.gca().set_xlim(x0,x1)
                x0=t0
                x1=t1
                #print("is_map?",obj.is_map)
                if obj.is_map:
                    if not style:
                        style="."
                    adata=self.iterplot(names,curves,varlist=var_list,style=style)
                else:
                    if not style:
                        style="-"
                    adata=self.dif_eq_plot(names,curves,varlist=var_list,v_slice=v_slice,polar=polar,style=style,scale=scale)
                    

                if anim:
                    for data in adata:
                        parinds.append(data)
                continue
            elif isinstance(obj,Inequality):# and "∈" in obj:
                obj=obj.eval(**at)
                if isinstance(obj,Inequality):
                    raise NotImplementedError("Inequalities in multiple variables haven't been implemented yet.")
                bounds,closed=self.extract_bounds(obj)
                label="Solution"
                for l,r in bounds:
                    orl,orr=l,r
                    if l=="empty":
                        continue
                    if l=="-∞":
                        if r=="∞":
                            l,r=x0,x1
                        else:
                            l=x0
                            if r<l:
                                l=r-1
                        #ax.scatter(l,0,c="none",edgecolors="b")
                    if r=="∞":
                        r=x1
                        if r<l:
                            d=l-r
                            r=l+1
                        #ax.scatter(r,0,c="none",edgecolors="b")
                    print(x0,x1,l,r)
                    x0=min(x0,l)
                    x1=max(x1,r)
                    print(x0,x1)
                    
                    ax.hlines(0,l,r,color="b",zorder=1000,label=label)
                    label=None
                    #ax.axhspan(l,r,color="b")
                    
                    if orl!="-∞":
                        ax.scatter(l,0,c="b" if closed else "none" ,edgecolors="b")
                    if orr!="∞":
                        ax.scatter(r,0,c="b" if closed else "none" ,edgecolors="b")
                continue
            elif isinstance(obj,Equation):
                if waveanim:
                    raise NotImplementedError("Animations for equations are not implemented yet")
                if ind in iinds:
                    vals=imagvals
                else:
                    for point in grid:
                        try:
                            newval=obj.eval(format=False,fractions=False,**{**{var:point},**at})
                            if isinstance(newval,str):
                                newval=None
                            if ind in rinds:
                                try:
                        
                                        imagvals.append(newval.imag)
                                        newval=newval.real
                                      
                                except AttributeError:
                                    imagvals.append(0)    
                                    pass
                        
                            else:
                                if isinstance(newval,Vector):
                                    newval=Vector(*(None if (isinstance(comp,Complex) and not comp.isreal()) else comp.real for comp in newval))
                        except ZeroDivisionError:
                            if isinstance(newval,Vector):
                                newval=Vector(*(None for comp in newval))
                            else:
                                newval=None
                        vals.append(newval)    
            else:
                if ind in iinds:
                    if waveanim:
                        wave_names.append(names[ind])
                        wavearr.append(imagtarr)
                    else:
                        vals=imagvals
                else:
                    if waveanim:
                        wave_names.append(names[ind])
                        if obj.numpiable():
                            tarr=np.array([np.ones(np.shape(grid))*obj.numcompose(**{**{spacevar:grid,timevar:t},**at}) for t in time_range])
                            if ind in rinds:
                                imagtarr=np.array(tarr.imag)
                                wavearr.append(tarr.real)
                            else:
                                wavearr.append(tarr)
                            #print(np.shape(wavearr))
                        else:
                            tarr=np.zeros((tden,len(grid)))
                            imagtarr=np.zeros((tden,len(grid)))
                            for tind,t in enumerate(time_range):
                                for sind,point in enumerate(grid):
                                    try:
                                        newval=obj.eval(**{**{spacevar:point,timevar:t},**at})
                                        if ind in rinds:
                                            try:
                                                imagtarr[tind,sind]=newval.imag
                                                newval=newval.real             
                                            except AttributeError:    
                                                imagvals.append(0)
                                    except ZeroDivisionError:
                                        newval=None
                                    tarr[tind,sind]=newval
                    else:
                        if obj.numpiable():
                            
                            vals=np.ones(np.shape(grid))*obj.numcompose(**{**{var:grid},**at})
                            
                            if ind in rinds:
                                imagvals=np.array(vals.imag)
                                vals=vals.real
                        else:
                            for point in grid:
                                try:
                                    newval=obj.eval(**{**{var:point},**at})
                                    if ind in rinds:
                                        try:
                                            imagvals.append(newval.imag)
                                            newval=newval.real
                                        except AttributeError:    
                                            imagvals.append(0)
                                        
                                    
                                except  ZeroDivisionError:
                                    newval=None
                                    if ind in rinds:
                                        imagvals.append(None)
                                vals.append(newval)
            #leg=re.sub("\^\((-?\d+)\)",lambda s: f"${{}}^{{{s[1]}}}$",names[ind])
            #leg=re.sub("\^\s*(-?\d+\s*)",lambda s: f"${{}}^{{{s[1]}}}$",leg)
            #leg="$x^x$"
            #leg="$x^{x+1}$"
            leg=names[ind]#"$"+self.nicepower(names[ind],lsymb="^{",rsymb="}")+"$"
            #print(leg)
            if isinstance(obj,Vector): #parametric plot
                #print(np.shape(vals))
                if not obj.get_vars() and obj.dim==2:
                    ax.plot([obj[1]],[obj[2]],".") #label=str(ind+1)+": "+leg+f", {var} $\in [ {x0:4g} , {x1:4g} ]$")
                    #xvals,yvals=zip(*vals)
                    #raise ValueError(str(xvals))
                    xvals,yvals=vals
                    #raise ValueError(str(vals))
                    self.plotdata.append([[obj[1]],[obj[2]]])
                    parinds.append(("I am a point",xvals,yvals))

                else:
                    if obj.numpiable():
                        vals=np.swapaxes(vals,0,1)
                    if polar:
                        yvals,xvals=zip(*vals)    
                    else:
                        xvals,yvals=zip(*vals)
                    if scale:
                        xvals=np.array(xvals);yvals=np.array(yvals)
                        xvals=xvals/max(abs(xvals))
                        yvals=yvals/max(abs(yvals))
                    if anim is None:
                        anim=True
                    if not waveanim:
                        if anim:
                            parinds.append((names[ind],xvals,yvals))
                        self.plotdata.append([xvals,yvals])
                        ax.plot(xvals,yvals,style,label=str(ind+1)+": "+leg+f", {var} $\in [ {x0:4g} , {x1:4g} ]$")
                    if len(varset)==1:
                        if var.strip() not in ("x","y"):
                            ax.set_xlabel("x")  
                            ax.set_ylabel("y")  
                        elif var.strip()=="x":
                            ax.set_xlabel("y")  
                            ax.set_ylabel("z")  
                        elif var.strip()=="y":
                            ax.set_xlabel("x")  
                            ax.set_ylabel("z")  
            else:
                if not waveanim:
                    if scale:
                        vals=np.array(vals)
                        vals=vals/max(abs(vals))
                    if anim:
                        if isinstance(obj,Equation):
                            dim=max(val.dim if isinstance(val,Vector) else 1 for val in vals)
                            
                            for i in range(1,dim+1):
                                parinds.append((names[ind],grid,[val if not isinstance(val,Vector) else (val[i] if val.dim>=i else None) for val in vals]))
                                
                        else:
                            parinds.append((names[ind],grid,vals))    
                    #print(len(grid),len(vals))
                    
                    self.plotdata.append([grid,vals])
                    ax.plot(grid,vals,style,label=str(ind+1)+": "+leg)
               # print("varset:",varset)
                ax.set_xlabel(", ".join(str(var) for var in sorted(varset)))
                last=chr(ord(sorted(varset)[-1][0]))
                
                if last=="z":
                    if "f" not in varset:
                        ax.set_ylabel("f")
                    elif "♥" not in varset:
                        ax.set_ylabel("♥")
                    else:
                        first=chr(ord(sorted(varset)[0][0])-1)
                        ax.set_ylabel(first)
                        
                else:
                    ax.set_ylabel(chr(ord(sorted(varset)[-1][0])+1))
        if is_open:
            for hline in hlines:
                ax.axhline(hline,lw=1,color="black")
            for vline in vlines:
                ax.axvline(vline,lw=1,color="black")
            if y0 is not None:
                plt.gca().set_ylim(bottom=y0)
            if y1 is not None:
                plt.gca().set_ylim(top=y1)
            #print(parinds)
            if not waveanim:
                ax=plt.gca() if hold else ax
                h,l=ax.get_legend_handles_labels()
                if len(l)>10:
                    h=h[:2]+[h[3],h[-3]]+h[-2:]
                    l=l[:2]+["."*len(l[1])]+["."*len(l[1])]+l[-2:]
                if legend:
                    fig.legend(h,l) #Zobrazíme maximálně 10 legend
            if anim:
                #print(parinds)
                if waveanim:
                    self.wave_anim(wave_names,grid,wavearr,t0,t1,fig,ax,polar=polar,style=style,equal=equal,save=save,
                    filename=filename,realtime=realtime,legend=legend,scale=scale,tlabel=timevar)
                else:
                    self.par_anim(parinds,x0,x1,polar=polar,style=style,equal=equal,pointanim=pointanim,save=save,
                    filename=filename,xlabel=plt.gca().get_xlabel(),ylabel=plt.gca().get_ylabel(),realtime=realtime,legend=legend,scale=scale)
            if equal:
                plt.gca().set_aspect("equal")
            if use_grid is None:
                use_grid=False
            if use_grid:
                plt.gca().grid(True)
            #animace:
            if loglog:
                plt.xscale("log")
                plt.yscale("log")
            if save:
                for fig_ind in plt.get_fignums():
                    plt.figure(fig_ind)
                    plt.savefig(filename+"1D"+str(fig_ind)+".png")
                    plt.close()
                #plt.savefig(filename+".png")
                
                matplotlib.use(self.backend)
            else:
              #  if self.backend=="Qt5Agg":
              #      plt.ion()
                #plt.ion()
                plt.show()      
                pass
        elif show2D:
            plt.show()
    
    def plot2(self,*objs,llim=None,ulim=None,pden=None,style="",save=False,filename=None,use_grid=None,contours=None,polar=None,equal=None,
    hold=False,anim=False,realtime=False,realimag=False,t0=0,t1=2*math.pi,tden=None,v_slice=None,legend=True,animt=None,**at):
        self.ac2=self.AnimControl()
        if contours==True:
            contours=20
        if save:
            matplotlib.use("Agg")
            if filename==None:
                filename="SavedResults/2DPlot"
            else:
                filename="SavedResults/"+filename
        else:
            #matplotlib.use(self.backend)
            matplotlib.use(self.backend)
        
        try:
            x0,y0=llim
        except TypeError:
            x0=y0=llim
        try:
            x1,y1=ulim
        except TypeError:
            x1=y1=ulim
        
        if style:
            cmap=style
            angcmap=style
        else: 
            cmap="plasma"
            angcmap="plasma" if realimag else "hsv"
        newobjs=[]
        
        for obj in objs:
            
            if isinstance(obj,Vector) and obj.dim>2:    
                if v_slice not in (None,"all"):
                    newobjs.append(Vector(obj[v_slice[0]],obj[v_slice[1]]))
                else:
                    for comb in it.combinations(obj,2):
                        newobjs.append(Vector(*comb))
            else:
                newobjs.append(obj)
        objs=newobjs
        print("|".join(str(obj) for obj in objs))
        for obj in objs:    
            if anim:
                animdata=[]
                
            if isinstance(obj,GeneralObject):
                obj=obj.eval(**at)
            if not isinstance(obj,GeneralObject):
                #obj=Monom(obj)
                obj=gen_obj.Const(gen_obj.Variable("x"),obj)

            vars=obj.get_vars()
            lv=len(vars)
            if lv==0:
                var1,var2="x","y"
                if anim:
                    t="t"
            elif lv==1:
                
                t=next(iter(vars))
                #print(t)
                var1,var2= (t,"y") if t!="y" else ("x","y")
                #print("polar...",polar)
                if polar is None:
                    if t=="r":
                        polar=True
                    elif t in ("f","fi","Φ","α","β","γ","φ","θ"):
                        polar=True
                        var1,var2=var2,var1
                elif polar is True:
                    if t in ("f","fi","Φ","α","β","γ","φ","θ"):
                        var1,var2=var2,var1
                if anim:
                    if "t" in (var1,var2):
                        if "s" in (var1,var2):
                            t="p"
                        else:
                            t="s"
                    else:
                        t="t"
                    
                #var1,var2=  "y","z"
                
            elif lv==2:
                ivars=iter(sorted(vars))
                var1,var2=ivars
                if polar is None:
                    if var2=="x" or var1=="y":
                        var1,var2=var2,var1
                    if any(symb in vars for symb in ("r","f","fi","Φ","α","β","γ","φ","θ")):
                        polar=True
                if polar is True:
                    if var2=="r" or var1 in ("f","fi","Φ","α","β","γ","φ","θ"):
                        var1,var2=var2,var1
                if anim:
                    if animt in vars:
                        t=animt
                        var1=next(iter(vars.difference({animt,})))
                     #   while var2 in {t,var2}:
                     #      var2=chr(ord(var2)+1) 
                    elif "t" in vars:
                        t="t"
                        var1=next(iter(vars.difference({"t",})))
                    else:
                        
                        if var1=="x":
                            t=var2
                        else:
                            t,var1=var1,var2
                    var2="_"#chr(min(ord(var1),ord(t))+1)
                    while var2 in {t,var1}:
                        var2=chr(ord(var2)+1)                     
                    if var1=="y":
                        var1,var2=var2,var1
                    #var2,t=chr(ord(sorted((var1,var2))[-1][0])+1),var2
                    

            else:
                if anim is False or lv>3:
                    raise NotImplementedError("Can only plot functions of at most two variables")
                else:
                    anim=True
                    if not animt in vars:
                        
                        t,var1,var2=sorted(vars) 
                    else:
                        t=animt
                        var1,var2=sorted(vars.difference({animt,}))

                    
            if anim:
                if t0 is None:
                    t0=0
                if t1 is None:
                    t1=2*np.pi
                if t1<t0:
                    t0,t1=t1,t0
                if tden is None or tden<0:
                    tden=100
                time_range=np.linspace(t0,t1,tden)
                #print(time_range)       

            if polar:
                if x0 is None:
                    x0=0
                if x1 is None:
                    x1=1
                if y0 is None:
                    y0=0
                if y1 is None:
                    y1=2*math.pi
            else:
                if x0 is None:
                    x0=-1
                if x1 is None:
                    x1=1
                if y0 is None:
                    y0=-1
                if y1 is None:
                    y1=1
            if x1<x0:
                x0,x1=x1,x0
            if y1<y0:
                y0,y1=y1,y0

            
            iscomplex=False
            isvector=False
            #if isinstance(obj,Complex):
            if obj.contains_complex():
                    iscomplex=True
                    if contours==None:
                        contours=20#True
            elif isinstance(obj,Vector):
                if obj.dim>3:
                    raise NotImplementedError("Can't plot vector fields in more than two dimensions")
                try:
                    obj1,obj2=(comp if isinstance(comp,GeneralObject) else Monom(comp) for comp in obj)
                except ValueError:
                    obj1=next(iter(obj))
                    obj1,obj2=obj1 if isinstance(obj1,GeneralObject) else Monom(obj1),Monom(0)
                isvector=True
                if contours==None:
                    contours=0#False
            if contours==None:
                contours=20#True

            #leg=re.sub("\^\((-?\d+)\)",lambda s: f"${{}}^{{{s[1]}}}$",str(obj))
            #leg=re.sub("\^\s*(-?\d+\s*)",lambda s: f"${{}}^{{{s[1]}}}$",leg)
            leg="$"+self.nicepower(str(obj),lsymb="^{",rsymb="}")+"$"
            if pden is None:
                if obj.numpiable():
                    if anim:
                        pden=250 if iscomplex else 125
                    else:
                        pden=500 if iscomplex else 250
                else:
                    if anim:
                        pden=50
                    else:
                        pden=125 if iscomplex else 250

            xgrid,ygrid=np.meshgrid(np.linspace(x0,x1,pden+1),np.linspace(y0,y1,pden+1))
            vals=np.array([ [(xgrid[i,j],ygrid[i,j]) for i,_ in enumerate(xgrid)] for j,_ in enumerate(xgrid[0])] )
            

            if not isvector:
                if iscomplex:
                    if not isinstance(obj,Complex):
                        obj=Complex(obj)
                    norm,ang=np.empty(np.shape(xgrid)),np.empty(np.shape(xgrid))
                    if realimag:
                        r=gen_obj.Real(obj)
                        i=gen_obj.Imag(obj)
                        onorm=r if isinstance(r,GeneralObject) else Monom(r)
                        oangle=i if isinstance(i,GeneralObject) else Monom(i)
                    else:
                        onorm=obj.norm() if isinstance(obj.norm(),GeneralObject) else Monom(obj.norm())
                        oangle=obj.angle() if isinstance(obj.angle(),GeneralObject) else Monom(obj.angle())
                    

                    if anim:
                        if onorm.numpiable() and oangle.numpiable():
                            #print("Numpiable!")

                            normanimdata=[onorm.numcompose(**{**{var1:xgrid,var2:ygrid,t:tind},**at}) for tind in time_range]
                            anganimdata=[oangle.numcompose(**{**{var1:xgrid,var2:ygrid,t:tind},**at}) for tind in time_range]
                            
                        else:
                            normanimdata=[]
                            anganimdata=[]
                            for tind in time_range:
                                norm=np.empty(np.shape(xgrid))
                                ang=np.empty(np.shape(xgrid))
                                for i,_ in enumerate(vals):
                                    for j,_ in enumerate(vals[:,0]):
                                        try:
                                            norm[j,i]=onorm.eval(**{**{var1:vals[i,j][0],var2:vals[i,j][1],t:tind},**at})
                                            ang[j,i]=oangle.eval(**{**{var1:vals[i,j][0],var2:vals[i,j][1],t:tind},**at})
                                        except ZeroDivisionError:
                                            norm[j,i]=None
                                            ang[j,i]=None
                                normanimdata.append(norm)
                                anganimdata.append(ang)
                        print("passing to anim")
                        obj_str="$"+self.nicepower(str(obj),lsymb="^{",rsymb="}")+"$"
                        titles=(f"Real part of {obj_str}",f"Imag. part of {obj_str}") if realimag else (f"Absolute value of {obj_str}",f"Phase of {obj_str}")
                        self.anim2d(xgrid,ygrid,normanimdata,time_range[0],time_range[-1],realtime=realtime,ac=self.ac2,hold=True,
                        polar=polar,style=style,xlabel=var1,ylabel=var2,tlabel=t,title=titles[0])
                        self.anim2d(xgrid,ygrid,anganimdata,time_range[0],time_range[-1],realtime=realtime,ac=self.ac2,hold=True,
                        polar=polar,style=style,xlabel=var1,ylabel=var2,tlabel=t,title=titles[1])
                        #plt.show()
                    else:
                        if onorm.numpiable() and oangle.numpiable():
                            print("Numpiable!")
                            norm=onorm.numcompose(**{**{var1:xgrid,var2:ygrid},**at})
                            ang=oangle.numcompose(**{**{var1:xgrid,var2:ygrid},**at})
                        else:
                            for i,_ in enumerate(vals):
                                for j,_ in enumerate(vals[:,0]):
                                    try:
                                        norm[j,i]=onorm.eval(**{**{var1:vals[i,j][0],var2:vals[i,j][1]},**at})
                                        ang[j,i]=oangle.eval(**{**{var1:vals[i,j][0],var2:vals[i,j][1]},**at})
                                    except ZeroDivisionError:
                                        norm[j,i]=None
                                        ang[j,i]=None
                        if realimag:
                            pref=("Real part","Imag. part")
                        else:
                            pref=("Absolute value","Phase")
                        
                        for ind,(valgrid,colmap) in enumerate(((norm,cmap),(ang,angcmap))):
                            fig,ax=plt.subplots()
                            if ind==1:
                                backend = matplotlib.get_backend()
                                if self.backend=="Qt5Agg":
                                    mngr = plt.get_current_fig_manager()
                                    geom = mngr.window.geometry()
                                    x,y,dx,dy = geom.getRect()
                                    mngr.window.setGeometry(x+dx/2, y+dy/2, dx, dy)
                            obj_str=self.nicepower(str(obj),lsymb="^{",rsymb="}")
                            fig.suptitle("{3} of $f({0},{1}) ={2}$".format(var1,var2,obj_str,pref[ind]))
                            if polar:
                                plt.axes(projection='polar')
                                im=plt.pcolormesh(ygrid,xgrid,valgrid,cmap=colmap)
                                cbar=plt.colorbar(im)
                                if contours:
                                    cset=np.linspace(np.min(valgrid),np.max(valgrid),contours)
                                    try:
                                        cnt=plt.contour(valgrid.transpose(),cset,linewidths=1,extent=(y0,y1,x0,x1),colors="black")
                                        cbar.add_lines(cnt)
                                    except ValueError as e:
                                        if str(e)=="Contour levels must be increasing":
                                            pass
                                        else:
                                            raise
                            else:
                                im=plt.pcolormesh(xgrid,ygrid,valgrid,cmap=colmap)
                                cbar=plt.colorbar(im)
                                if contours:
                                    cset=np.linspace(np.min(valgrid),np.max(valgrid),contours)
                                    try:
                                        cnt=plt.contour(valgrid,cset,linewidths=1,extent=(x0,x1,y0,y1),colors="black")
                                        cbar.add_lines(cnt)
                                    except ValueError as e:
                                        if str(e)=="Contour levels must be increasing":
                                            pass
                                        else:
                                            raise
                            if equal!=False:
                                ax.set_aspect("equal")
                            if save:
                                plt.savefig(filename+"_"+pref[ind]+str(self.plot_ind)+".png")
                                self.plot_ind+=1
                                plt.close()
                                #plt.close(fig)
                                matplotlib.use(self.backend)
                    
                else:
                    if anim:
                        if obj.numpiable():
                            print("Numpiable!")
                            animdata=[obj.numcompose(**{**{var1:xgrid,var2:ygrid,t:tind},**at}) for tind in time_range]
                            
                        else:
                            for tind in time_range:
                                valgrid=np.empty(np.shape(xgrid))
                                for i,_ in enumerate(vals):
                                    for j,_ in enumerate(vals[:,0]):
                                        try:
                                            valgrid[j,i]=obj.eval(**{**{var1:vals[i,j][0],var2:vals[i,j][1],t:tind},**at})
                                            
                                        except ZeroDivisionError:
                                            #print("error",obj,valgrid[j,i])
                                            valgrid[j,i]=None    
                                animdata.append(valgrid)
                     #   print("hold",hold)
                        self.anim2d(xgrid,ygrid,animdata,time_range[0],time_range[-1],realtime=realtime,ac=self.ac2,
                        hold=True,polar=polar,style=style,xlabel=var1,ylabel=var2,tlabel=t,title="$"+self.nicepower(str(obj),lsymb="^{",rsymb="}")+"$")
                        #plt.show()

                    else:
                        fig,ax=plt.subplots()
                        if obj.numpiable():
                            valgrid=obj.numcompose(**{**{var1:xgrid,var2:ygrid},**at})
                        else:
                            print("Not Numpiable!")
                            valgrid=np.empty(np.shape(xgrid))
                            
                            for i,_ in enumerate(vals):
                                for j,_ in enumerate(vals[:,0]):
                                    try:
                                        valgrid[j,i]=obj.eval(**{**{var1:vals[i,j][0],var2:vals[i,j][1]},**at})
                                        #print(obj,valgrid[j,i])
                                        #input()
                                    except ZeroDivisionError:
                                        #print("error",obj,valgrid[j,i])
                                        valgrid[j,i]=None
                    
                        if polar:
                            plt.axes(projection='polar')
                            im=plt.pcolormesh(ygrid,xgrid,valgrid,cmap=cmap)
                            cbar=plt.colorbar(im)
                            if contours:
                                cset=np.linspace(np.min(valgrid),np.max(valgrid),contours)
                                cnt=plt.contour(valgrid.transpose(),cset,linewidths=1,extent=(y0,y1,x0,x1),colors="black")
                                cbar.add_lines(cnt)
                        else:
                            im=plt.pcolormesh(xgrid,ygrid,valgrid,cmap=cmap)
                            cbar=plt.colorbar(im)
                            if contours:
                                cset=np.linspace(np.min(valgrid),np.max(valgrid),contours)
                                cnt=plt.contour(valgrid,cset,linewidths=1,extent=(x0,x1,y0,y1),colors="black")
                                cbar.add_lines(cnt)
                        plt.title("Plot of f({0},{1}) ={2}".format(var1,var2,leg))
            else:
               
                
                if anim:
                    if obj.numpiable():
                        #print("Numpiable!")
                        animdata=[obj.numcompose(**{**{var1:xgrid,var2:ygrid,t:tind},**at}) for tind in time_range]
                        #print(animdata)
                    else:
                        animdata=[]
                        for tind in time_range:
                            valgrid=np.empty(np.shape(xgrid))
                            for i,_ in enumerate(vals):
                                for j,_ in enumerate(vals[:,0]):
                                    try:
                                        valgrid[j,i]=obj.eval(**{**{var1:vals[i,j][0],var2:vals[i,j][1],t:tind},**at})
                                        
                                    except ZeroDivisionError:
                                        #print("error",obj,valgrid[j,i])
                                        valgrid[j,i]=None    
                            animdata.append(valgrid)
                    
                    self.anim2d(xgrid,ygrid,np.array(animdata),time_range[0],time_range[-1],realtime=realtime,ac=self.ac2,
                    hold=True,polar=polar,vecfield=True,style=style,xlabel=var1,ylabel=var2,tlabel=t,title=str(obj))

                else:
                    fig,ax=plt.subplots()
                    if obj.numpiable():
                        #print("Numpiable!")
                        valgrids=obj.numcompose(**{**{var1:xgrid,var2:ygrid},**at})
                    else:
                        valgrids=[np.empty(np.shape(xgrid)),np.empty(np.shape(xgrid))]
                        for i,_ in enumerate(vals):
                            for j,_ in enumerate(vals[:,0]):
                                try:
                                    valgrids[0][j,i]=obj1.eval(**{**{var1:vals[i,j][0],var2:vals[i,j][1]},**at})
                                    valgrids[1][j,i]=obj2.eval(**{**{var1:vals[i,j][0],var2:vals[i,j][1]},**at})
                                except ZeroDivisionError:
                                    valgrid[j,i]=None
                    #print(var1,var2,np.max(xgrid),np.max(ygrid),xgrid[0]])        
                    if polar:
                        plt.axes(projection='polar')
                        #valgrids[0],valgrids[1]=valgrids[1],valgrids[0]#*ygrid
                        norm=np.sqrt(valgrids[0]**2+valgrids[1]**2)
                        im=plt.pcolormesh(ygrid,xgrid,norm,cmap="plasma")
                        
                        lx,ly=np.shape(ygrid)
                        decim=30
                        dy,dx=ly//decim,lx//decim
                        dxgrid,dygrid=[xgrid[::dx,::dy],ygrid[::dx,::dy]]
                        dvalgrids=[valgrids[0][::dx,::dy],valgrids[1][::dx,::dy]]
                        plt.quiver(dygrid,dxgrid,dvalgrids[0]*np.cos(dygrid)-dvalgrids[1]*np.sin(dygrid), dvalgrids[0]*np.sin(dygrid)+dvalgrids[1]*np.cos(dygrid))
                        if contours:
                            cset=np.linspace(np.min(norm),np.max(norm),contours)
                            cnt=plt.contour(norm.transpose(),cset,extent=(y0,y1,x0,x1))
                    else:
                        norm=np.sqrt(valgrids[0]**2+valgrids[1]**2)
                        

                        im=plt.pcolormesh(xgrid,ygrid,norm,cmap=cmap)
                        cbar=plt.colorbar(im)
                        lx,ly=np.shape(ygrid)
                        decim=30
                        dy,dx=ly//decim,lx//decim
                        plt.quiver(xgrid[::dx,::dy],ygrid[::dx,::dy],valgrids[0][::dx,::dy],valgrids[1][::dx,::dy])
                        if contours:
                            cset=np.linspace(np.min(norm),np.max(norm),contours)
                            cnt=plt.contour(norm,cset,extent=(x0,x1,y0,y1))
            
            
                
                    plt.title("Plot of f({0},{1}) ={2}".format(var1,var2,leg))
                
            if not anim:
                plt.gca().set_xlabel(var1)
                plt.gca().set_ylabel(var2)
                if equal!=False:
                    plt.gca().set_aspect("equal")
                if use_grid is None:
                    use_grid=False
                if use_grid:
                    plt.gca().grid(True)
        
        if save:
            plt.savefig(filename+str(self.plot_ind)+".png")
            plt.close()
            self.plot_ind+=1
            matplotlib.use(self.backend)
        else:
            if self.backend=="Qt5Agg":
                if not anim:
                    plt.ion()
                else:
                    plt.ioff()
            
            if not hold:
                plt.show()      
    
    class AnimControl: #Stavové proměnné animace/animací  

        def __init__(self):
            self.pause=False
            self.direction=1
            self.frame=self.ref_frame=self.ref_time=0

    def par_anim(self,data,t0,t1,polar=False,style="",equal=False,save=False,pointanim=False,filename=None,realtime=False,legend=True,xlabel="x",ylabel="y",scale=False):
        
        from matplotlib.animation import FuncAnimation
        if save:
            if filename==None:
                filename="SavedResults/Anim"
            else:
                filename="SavedResults/"+filename    
        if not data:
            return
        
        if polar:
            fig=plt.figure()
            ax = fig.add_subplot(111, projection='polar')
        else:
            fig,ax=plt.subplots()

        if  matplotlib.get_backend()=="Qt5Agg":
            mngr = plt.get_current_fig_manager()
            geom = mngr.window.geometry()
            x,y,dx,dy = geom.getRect()
            mngr.window.setGeometry(x+dx/2, y+dy/2, dx, dy)
        
        names,xdata,ydata=zip(*data)
        #print(names,xdata)
        xdata=[tuple(val.real if isinstance(val,(complex,Complex)) else val for val in tup) for tup in xdata]
        ydata=[tuple(val.real if isinstance(val,(complex,Complex)) else val for val in tup) for tup in ydata]


        duration=t1-t0
        len_t=len(xdata[0])
        dt=duration/len_t
        if pointanim:
            #lines=[ax.plot([],[],style,marker="o",label=name)[0] for name in names]    
            lines=[ax.plot([],[],".")[0] if name=="I am a point" else ax.plot([],[],style,label=name,marker="o")[0]  for name in names]
        else:
            lines=[ax.plot([],[],".")[0] if name=="I am a point" else ax.plot([],[],style,label=name)[0]  for name in names]
        
        if polar:
             x_min,x_max=0,2*math.pi
        else:
            x_min=np.min(xdata)
            x_max=np.max(xdata)
        y_min=np.min(ydata)
        y_max=np.max(ydata)
        time_text = ax.text(.7, .8, '', fontsize=15,transform=ax.transAxes)
        
        def init():
            #ax.clear()
            for line in lines:
                line.set_data([],[])
            ax.set_xlim(x_min,x_max)
            ax.set_ylim(y_min,y_max)
            time_text.set_text("")
            self.ac.frame=0
            self.ac.ref_time=time.time()+t0
            h,l=plt.gca().get_legend_handles_labels()
            if len(l)>10:
                h=h[:2]+[h[3],h[-3]]+h[-2:]
                l=l[:2]+["."*len(l[1])]+["."*len(l[1])]+l[-2:]
            if legend:
                fig.legend(h,l) #Zobrazíme maximálně 10 legend
            ax.set_xlabel(xlabel);ax.set_ylabel(ylabel)
            return lines+[time_text]
        
        def update(_):
            if not self.ac.pause:
                if realtime:
                    self.ac.frame=(self.ac.ref_frame+self.ac.direction*int((time.time()-self.ac.ref_time)/(dt))) %len_t # k uloženému číslu framu přidá číslo odpovídající času uplynulému od posledního uložení
                else:
                    self.ac.frame=(self.ac.frame+self.ac.direction)%len_t
                    
            #if self.ac.frame==len_t:
            #    self.ac.frame=0
            #elif self.ac.frame==-1:
            #    self.ac.frame=len_t-1
            for lind,line in enumerate(lines):    
                line.set_data(xdata[lind][:self.ac.frame],ydata[lind][:self.ac.frame])
            time_text.set_text("{tlabel} = {val:.3f}".format(tlabel=xlabel,val=t0+(self.ac.frame)*dt))
            return lines+[time_text]

        def pointupdate(_):
            if not self.ac.pause:
                if realtime:
                    self.ac.frame=(self.ac.ref_frame+self.ac.direction*int((time.time()-self.ac.ref_time)/(dt))) %len_t
                else:
                    self.ac.frame=(self.ac.frame+self.ac.direction)%len_t
            #if self.ac.frame==len_t:
            #    self.ac.frame=0
            #elif self.ac.frame==-1:
            #    self.ac.frame=len_t-1

            for lind,line in enumerate(lines):    
                line.set_data([xdata[lind][self.ac.frame]],[ydata[lind][self.ac.frame]])
            time_text.set_text("{tlabel} = {val:.3f}".format(tlabel=xlabel,val=t0+(self.ac.frame)*dt))
            return lines+[time_text]
        
        def onClick(event):
            if event.button==1:
                if not self.ac.pause:
                    self.ac.ref_frame=self.ac.frame
                else:
                    self.ac.ref_time=time.time()+t0
                self.ac.pause = not self.ac.pause
            
            elif event.button==3:
                self.ac.ref_time=time.time()+t0
                self.ac.ref_frame=self.ac.frame
                self.ac.direction=-self.ac.direction

        def onPress(event):
            
            if event.key in [" "]:
                if not self.ac.pause:
                    self.ac.ref_frame=self.ac.frame
                else:
                    self.ac.ref_time=time.time()+t0
                self.ac.pause = not self.ac.pause
            elif event.key=="right":
                if self.ac.direction>0:
                    self.ac.direction=2*self.ac.direction
                else:
                    if self.ac.direction<-1:
                        self.ac.direction=min(-1,self.ac.direction//2)
                    elif self.ac.direction==-1:
                        self.ac.direction=0
                    else:
                        self.ac.direction=1
            elif event.key=="left":
                if self.ac.direction<0:
                    self.ac.direction=2*self.ac.direction
                else:
                    if self.ac.direction>1:
                        self.ac.direction=max(1,self.ac.direction//2)
                    elif self.ac.direction==1:
                        self.ac.direction=0
                    else:
                        self.ac.direction=-1
            elif event.key=="up":
                self.ac.direction=1 if self.ac.direction>0 else -1
            elif event.key=="down":
                self.ac.frame=0
                self.ac.ref_frame=0
          
        def onScroll(event):
            self.ac.frame=(self.ac.frame+int(event.step))%len_t
            self.ac.ref_frame=self.ac.frame

        fig.canvas.mpl_connect('button_press_event', onClick)
        fig.canvas.mpl_connect('scroll_event', onScroll)
        fig.canvas.mpl_connect('key_press_event',onPress)

        if pointanim:
            self.anim=FuncAnimation(fig,pointupdate,init_func=init,blit=True,interval=1,repeat=True,save_count=len_t)
        else:
            self.anim=FuncAnimation(fig,update,init_func=init,blit=True,interval=1,repeat=True,save_count=len_t)
        if save:
            self.anim.save(filename+".mp4",fps=len_t//duration)
            """plt.rcParams['animation.ffmpeg_path'] = 'C:/FFmpeg/bin/ffmpeg.exe'
            FFwriter = matplotlib.animation.FFMpegWriter(fps=100)
            self.anim.save(filename+".mp4",writer=FFwriter)
            """
        if equal:
                plt.gca().set_aspect("equal")
                
    def wave_anim(self,names,grid,ydata,t0,t1,fig,ax,polar=False,style="",equal=False,save=False,filename=None,realtime=False,legend=True,scale=False,tlabel="t"):
        
        from matplotlib.animation import FuncAnimation
        
        if save:
            if filename==None:
                filename="SavedResults/Anim"
            else:
                filename="SavedResults/"+filename    
        """
        if polar:
            fig=plt.figure()
            ax = fig.add_subplot(111, projection='polar')
        else:
            fig,ax=plt.subplots()
        """
        if  matplotlib.get_backend()=="Qt5Agg":
            mngr = plt.get_current_fig_manager()
            geom = mngr.window.geometry()
            x,y,dx,dy = geom.getRect()
            mngr.window.setGeometry(int(x+dx//2), int(y+dy//2), int(dx), int(dy))
        
        #ydata=[tuple(val.real if isinstance(val,(complex,Complex)) else val for val in tup) for tup in ydata]
        

        duration=t1-t0
        len_t=len(ydata[0])
        dt=duration/len_t
        lines=[ax.plot([],[],style,lw=2,label=name)[0] for name in names]
        if scale:
            for i,dati in enumerate(ydata):
                yspacemax=np.max(abs(dati),axis=1)
            #print(np.shape(np.swapaxes(ydata,1,2)),np.shape(yspacemax))
                dati=np.swapaxes(dati,0,1)/yspacemax
                ydata[i]=np.swapaxes(dati,0,1)
        if polar:
             x_min,x_max=0,2*math.pi
        else:
            x_min=np.min(grid)
            x_max=np.max(grid)
        y_min=np.min(ydata)
        y_max=np.max(ydata)
        
        time_text = ax.text(.7, .8, '', fontsize=15,transform=ax.transAxes)
     #   print(np.shape(ydata),np.shape(names),duration,len_t,dt)
        def init():
            #ax.clear()
            for line in lines:
                line.set_data([],[])
            ax.set_xlim(x_min,x_max)
            ax.set_ylim(y_min,y_max)
            time_text.set_text("")
            self.ac.frame=0
            self.ac.ref_time=time.time()+t0
            h,l=plt.gca().get_legend_handles_labels()
            if len(l)>10:
                h=h[:2]+[h[3],h[-3]]+h[-2:]
                l=l[:2]+["."*len(l[1])]+["."*len(l[1])]+l[-2:]
            if legend:
                fig.legend(h,l) #Zobrazíme maximálně 10 legend
            return lines+[time_text]
        
        def update(_):
            if not self.ac.pause:
                if realtime:
                    self.ac.frame=(self.ac.ref_frame+self.ac.direction*int((time.time()-self.ac.ref_time)/(dt))) %len_t # k uloženému číslu framu přidá číslo odpovídající času uplynulému od posledního uložení
                else:
                    self.ac.frame=(self.ac.frame+self.ac.direction)%len_t
                    
          #  if self.ac.frame==len_t:
          #      self.ac.frame=0
          #  elif self.ac.frame==-1:
          #      self.ac.frame=len_t-1
            for lind,line in enumerate(lines):    
                line.set_data(grid,ydata[lind][self.ac.frame,:])
            time_text.set_text("{tlabel} = {val:.3f}".format(tlabel=tlabel,val=t0+(self.ac.frame)*dt))
            return lines+[time_text]

        def onClick(event):
            if event.button==1:
                if not self.ac.pause:
                    self.ac.ref_frame=self.ac.frame
                else:
                    self.ac.ref_time=time.time()+t0
                self.ac.pause = not self.ac.pause
            
            elif event.button==3:
                self.ac.ref_time=time.time()+t0
                self.ac.ref_frame=self.ac.frame
                self.ac.direction=-self.ac.direction

        def onPress(event):
            
            if event.key in [" "]:
                if not self.ac.pause:
                    self.ac.ref_frame=self.ac.frame
                else:
                    self.ac.ref_time=time.time()+t0
                self.ac.pause = not self.ac.pause
            elif event.key=="right":
                if self.ac.direction>0:
                    self.ac.direction=2*self.ac.direction
                else:
                    if self.ac.direction<-1:
                        self.ac.direction=min(-1,self.ac.direction//2)
                    elif self.ac.direction==-1:
                        self.ac.direction=0
                    else:
                        self.ac.direction=1
            elif event.key=="left":
                if self.ac.direction<0:
                    self.ac.direction=2*self.ac.direction
                else:
                    if self.ac.direction>1:
                        self.ac.direction=max(1,self.ac.direction//2)
                    elif self.ac.direction==1:
                        self.ac.direction=0
                    else:
                        self.ac.direction=-1
            elif event.key=="up":
                self.ac.direction=1 if self.ac.direction>0 else -1
            elif event.key=="down":
                self.ac.frame=0
                self.ac.ref_frame=0
        def onScroll(event):
            self.ac.frame=(self.ac.frame+int(event.step))%len_t
            if realtime:
                self.ac.ref_frame=self.ac.frame

        fig.canvas.mpl_connect('button_press_event', onClick)
        fig.canvas.mpl_connect('scroll_event', onScroll)
        fig.canvas.mpl_connect('key_press_event',onPress)

        
        self.anim2=FuncAnimation(fig,update,init_func=init,blit=True,interval=1,repeat=True,save_count=len_t)
        if save:
            self.anim2.save(filename+".mp4",fps=len_t//duration)
            """plt.rcParams['animation.ffmpeg_path'] = 'C:/FFmpeg/bin/ffmpeg.exe'
            FFwriter = matplotlib.animation.FFMpegWriter(fps=100)
            self.anim.save(filename+".mp4",writer=FFwriter)
            """
        if equal:
                plt.gca().set_aspect("equal")

    def anim2d(self,xgrid,ygrid,f,t0,t1,polar=False,style="",save=False,pointanim=False,filename=None,
    realtime=True,hold=False,ac=None,vecfield=False,xlabel="x",ylabel="y",tlabel="t",title=None): 
    # realtime=True => animace bude trvat (minimálně) t1-t0 s
    # hold=True => Animace se jen připraví, ale nespustí se
      #  print(np.min(f),np.max(f))
        from matplotlib.animation import FuncAnimation
        if save and realtime:
            realtime=False 
       # print(t0,t1)
        if ac is None: # Pokud bychom nechtěli synchronizovat
            ac=self.ac2
      #  print("Frame:",ac.frame,ac.ref_frame)
        if polar:
            fig=plt.figure()
            ax = fig.add_subplot(111, projection='polar')
        else:
            fig,ax=plt.subplots()
        if title is not None:
            #print(title)
            fig.suptitle(f"f({xlabel},{ylabel},{tlabel}) = {title}")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        #names,xdata,ydata=zip(*data)  # zip(*y) je inverze k zip(x), tj. místo 1 argumentu 'data' by se rovnou mohly vkládat 3 argumenty 'names','xdata','ydata'
        
        if t0>t1:   # t0 a t1,max_t,len_t,dt... jen kontrolují zobrazování a rychlost animace
            t0,t1=t1,t0
        if style:
            cmap=style
        else: 
            cmap="plasma"
        duration=t1-t0
        len_t=len(f)
        dt=duration/len_t
        if polar:
            xgrid,ygrid=ygrid,xgrid
        x_min=np.min(xgrid)
        x_max=np.max(xgrid)
        y_min=np.min(ygrid)
        y_max=np.max(ygrid)

        if vecfield:
            f=np.swapaxes(f,0,1)#f.reshape((2,len_t,*np.shape(xgrid)))
            f=f.reshape((2,len_t,*(np.shape(xgrid))))
            xcomp=f[0,:,:,:]
            ycomp=f[1,:,:,:]
            norm=np.sqrt(xcomp**2+ycomp**2)
            lx,ly=np.shape(ygrid)
            decim=30
            dy,dx=ly//decim,lx//decim
            dxgrid,dygrid=[xgrid[::dx,::dy],ygrid[::dx,::dy]]
            dvalgrids=f[:,:,::dx,::dy]
            f_min,f_max=np.min(norm),np.max(norm)
        else:
            f_min,f_max=np.min(f),np.max(f)
        global pc,qui
        #pc=ax.pcolormesh(np.empty(np.shape(xgrid)),np.empty(np.shape(xgrid)),np.empty(np.shape(xgrid)),label="2d plot",cmap=cmap)
        pc=ax.pcolormesh([[],[]],label="2d plot",cmap=cmap)
        if vecfield:
            qui=ax.quiver(np.empty(np.shape(dxgrid)),np.empty(np.shape(dygrid)),np.empty(np.shape(dvalgrids[0,ac.frame,:,:])), 
        np.empty(np.shape(dvalgrids[1,ac.frame,:,:])))
        time_text = ax.text(.7, .8, '', fontsize=15,transform=ax.transAxes) #transform= volí typ souřadnic, v tomhle případě relativně k plotovacímu rámečku 

        def init():
            global pc,qui
            ax.set_xlim(x_min,x_max)
            ax.set_ylim(y_min,y_max)
            time_text.set_text("")        
            
            ac.frame=0
            ac.ref_time=time.time()+t0
            if vecfield:
                qui=ax.quiver(np.empty(np.shape(dxgrid)),np.empty(np.shape(dygrid)),np.empty(np.shape(dxgrid)),np.empty(np.shape(dxgrid)))
                return [pc,qui,time_text]
            else:
                return [pc,time_text]
        
        def update(*_): # FuncAnim do update vkládá argument s číslem obrázku, já používám vlastní číslování (ac.frame), protože:
                    #  1) v případě více animací budou všechny synchronizované 2) tak můžu animaci ovládat přes myš 3) se mi se zapnutým blittingem jinak nepodařilo rozumně zprovoznit opakující se animace 
            global pc,qui
            if not ac.pause:
                if realtime:
                    ac.frame=(ac.ref_frame+ac.direction*int((time.time()-ac.ref_time)/(dt))) %len_t # k uloženému číslu framu přidá číslo odpovídající času uplynulému od posledního uložení
                else:
                    ac.frame=(ac.frame+ac.direction)%len_t
                    
          #  if ac.frame==len_t:
          #      ac.frame=0
          #  elif ac.frame==-1:
          #      ac.frame=len_t-1
            
            pc=ax.pcolormesh(xgrid,ygrid,f[ac.frame],label="2d plot",vmin=f_min,vmax=f_max,cmap=cmap)
            time_text.set_text("{0} = {1:.3f}".format(tlabel,t0+(ac.frame)*dt))
            return [pc,time_text]

        def update_vec_cart(*_):
            global pc,qui
            if not ac.pause:
                if realtime:
                    ac.frame=(ac.ref_frame+ac.direction*int((time.time()-ac.ref_time)/(dt))) %len_t # k uloženému číslu framu přidá číslo odpovídající času uplynulému od posledního uložení
                else:
                    ac.frame=(ac.frame+ac.direction)%len_t
         #   if ac.frame==len_t:
         #       ac.frame=0
         #   elif ac.frame==-1:
         #       ac.frame=len_t-1
            pc=ax.pcolormesh(xgrid,ygrid,norm[ac.frame,:,:],label="2d plot",vmin=0,vmax=f_max,cmap=cmap)
            #plt.quiver(dygrid,dxgrid,dvalgrids[0]*np.cos(dygrid)-dvalgrids[1]*np.sin(dygrid), dvalgrids[0]*np.sin(dygrid)+dvalgrids[1]*np.cos(dygrid))
            qui=ax.quiver(dxgrid,dygrid,dvalgrids[0,ac.frame,:,:], dvalgrids[1,ac.frame,:,:])
            time_text.set_text("{0} = {1:.3f}".format(tlabel,t0+(ac.frame)*dt))
            return [pc,qui,time_text]

        def update_vec_polar(*_):
            global pc,qui
            if not ac.pause:
                if realtime:
                    ac.frame=(ac.ref_frame+ac.direction*int((time.time()-ac.ref_time)/(dt))) %len_t # k uloženému číslu framu přidá číslo odpovídající času uplynulému od posledního uložení
                else:
                    ac.frame=(ac.frame+ac.direction)%len_t
         #   if ac.frame==len_t:
          #      ac.frame=0
          #  elif ac.frame==-1:
          #      ac.frame=len_t-1
            pc=ax.pcolormesh(xgrid,ygrid,norm[ac.frame,:,:],label="2d plot",vmin=0,vmax=f_max,cmap=cmap)
            qui=ax.quiver(dxgrid,dygrid,dvalgrids[0,ac.frame,:,:]*np.cos(dxgrid)-dvalgrids[1,ac.frame,:,:]*np.sin(dxgrid), dvalgrids[0,ac.frame,:,:]*np.sin(dxgrid)+dvalgrids[1,ac.frame,:,:]*np.cos(dxgrid))
            
            time_text.set_text("{0} = {1:.3f}".format(tlabel,t0+(ac.frame)*dt))
            return [pc,qui,time_text]

        def onClick(event):
            
            if event.button==1:
                if not ac.pause:
                    ac.ref_frame=ac.frame
                else:
                    ac.ref_time=time.time()+t0
                ac.pause = not ac.pause
                
            elif event.button==3:
                ac.ref_time=time.time()+t0
                ac.ref_frame=ac.frame
                ac.direction=-ac.direction

        def onPress(event):
            if event.key in [" "]:
                if not ac.pause:
                    ac.ref_frame=ac.frame
                else:
                    ac.ref_time=time.time()+t0
                ac.pause = not ac.pause
            elif event.key=="right":
                if ac.direction>0:
                    ac.direction=2*ac.direction
                else:
                    if ac.direction<-1:
                        ac.direction=min(-1,ac.direction//2)
                    elif ac.direction==-1:
                        ac.direction=0
                    else:
                        ac.direction=1
            elif event.key=="left":
                if ac.direction<0:
                    ac.direction=2*ac.direction
                else:
                    if ac.direction>1:
                        ac.direction=max(1,ac.direction//2)
                    elif ac.direction==1:
                        ac.direction=0
                    else:
                        ac.direction=-1
            elif event.key=="up":
                ac.direction=1 if ac.direction>0 else -1
            elif event.key=="down":
                ac.frame=0
                ac.ref_frame=0

        def onScroll(event):
            
            ac.frame=(ac.frame+int(event.step))%len_t
            if realtime:
                ac.ref_frame=ac.frame

        fig.canvas.mpl_connect('button_press_event', onClick)
        fig.canvas.mpl_connect('scroll_event', onScroll)
        fig.canvas.mpl_connect('key_press_event',onPress)

        if vecfield:
            if polar:
                self.anim=FuncAnimation(fig,func=update_vec_polar,init_func=init,blit=True,interval=0,repeat=True,save_count=len_t)
            else:
                self.anim=FuncAnimation(fig,func=update_vec_cart,init_func=init,blit=True,interval=0,repeat=True,save_count=len_t)
        else:
            self.animlist.append(FuncAnimation(fig,func=update,init_func=init,blit=True,interval=0,repeat=True,save_count=len_t))
            #self.anim=(FuncAnimation(fig,func=update,init_func=init,blit=True,interval=0,repeat=True,save_count=len_t))
            #self.animlist.append(self.anim)

        if save:
            if filename is None:
                filename="SavedResults/Anim"
            else:
                filename="SavedResults/"+filename     
            self.anim.save(filename+".mp4",fps=len_t//duration)
        else:
            
            if not hold:
                plt.show() 
            else:
                plt.show(block=False)   

    def iterplot(self,names,pointset,varlist=["x","y"],realtime=False,style="."):
        if not varlist:
            xlabel,ylabel="x","y"
        else:
            if len(varlist)<2:
                xlabel=varlist[0]
                ylabel="y" if xlabel!="y" else "z"
            else:
                xlabel,ylabel=varlist
        animdata=[]
        for name,points in zip(names,pointset):
            first=next(iter(points))
            dim=first.dim
            
            if dim==1:
                ax=plt.gca()
                if points.contains_complex():
                    points=[point[1] for point in points]
                    #print(points)
                    
                    xdata,ydata=zip(*((c.real,c.imag) if isinstance(c,Complex) else (c,0) for c in points))
                    
                    xlabel,ylabel=f"Real {xlabel}",f"Imag {xlabel}"
                else:
                    points=[point[1] for point in points]
                    xdata,ydata=list(range(len(points))),points
                
                #if anim:
                #    self.par_anim([("Map",list(range(points.dim)),points)],t0=0,t1=points.dim/100,style=".",realtime=realtime)
                #else:
                
            elif dim==2:
                ax=plt.gca()
                xdata,ydata=zip(*points)
            self.plotdata.append([xdata,ydata])
            ax.plot(xdata,ydata,style,label=name)
            animdata.append([name,xdata,ydata])
            ax.set_xlabel(xlabel);ax.set_ylabel(ylabel)
                #if anim:
                #    self.par_anim([("Map",*zip(*points))],t0=0,t1=points.dim/100,style=".",realtime=realtime)
        #plt.show()
        return animdata

    def dif_eq_plot(self,names,curves,*args,varlist=None,tlabel="t",v_slice=None,polar=False,style=".",scale=False,**kwargs):
        ax=plt.gca()#subplots()
        animdata=[]
        #print(names,curves)
        for name,curve in zip(names,curves):
            if isinstance(curve,np.ndarray):
                if v_slice is None:
                    _,m=np.shape(curve)
                    
                    if m==2:
                        v_slice=[0,1]
                        
                    else:
                        v_slice=[1,2]
                        
                
                xdata,ydata=curve[:,v_slice[0]],curve[:,v_slice[1]]
                
                if polar:
                    xdata,ydata=ydata,xdata
                if scale:
                    ydata=ydata/max(np.abs(ydata))
                self.plotdata.append([xdata,ydata])
                ax.plot(xdata,ydata,style,label=name)
                animdata.append([name,xdata,ydata])
        if varlist is None:
            xlabel,ylabel="x","y"
        else:
            #print(v_slice,varlist)
            #print(varlist,v_slice)
            xlabel,ylabel=varlist[v_slice[0]],varlist[v_slice[1]]
        ax.set_xlabel(xlabel);ax.set_ylabel(ylabel)
        #ax.legend()
        #if anim:
        #    self.par_anim(animdata,0,10,realtime=True)
        return animdata
#        plt.show()



