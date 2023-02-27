from rdkit import Chem
from rdkit.Chem import AllChem
import os
import math
from openbabel import pybel
import random
import numpy


def once():
    global core,layer,complexity
    complexity=0 
    do_cyclization(core)
    return not core_grp_fill()
    
def core_grp_fill():
    global unique_grp,core,layer,grp_choice   
    
    for grp in unique_grp:
        symm=symm_cast(grp)
        grp_choice=[]
        for i in range(len(symm)):
            k=symm.index(symm[i])
            if(k==i):
                flag=fill_grp(grp,-1)
            else:
                flag=fill_grp(grp,k)
            if(flag==1):
                return 1
            if(flag==2):
                break
    return 0

def symm_cast(s):
    global grp_symm,grp_cast
    n=len(grp_cast)
    symm=[ grp_symm[i] for i in range(n) if grp_cast[i]==s]
    return symm

def fill_grp(s,k):
    global core_type,layer,complexity,grp_choice
    u=layer.find("("+s,0)
    if(u<0):
        return 2
    else:
        v=layer.find(")",u)
        m=layer.find("%",u+1,v)
        if(m>=0):
             k=-1
        if(k<0):
            q=select_core_grp(s)
        else:
            q=grp_choice[k]
        if(q!="[H]"):
            complexity+=1
        q=type_fill(q)
        grp_choice.append(q)
        if(m>=0):
            z=layer[m:v]
            q=label_cyclic_type(q,z)
        else:
            if "*" in q:
                return 1
        if(len(q)):
            layer=layer[0:u+1]+q+layer[v:] 
            return 0
        else:
            return 1

def type_fill(s) :
    global type_choice
    type_choice=s
    while substitute_type_once(0)+substitute_type_once(1):
        continue   
    return type_choice     

def substitute_type_once(a):
    global type_choice,grp_type,wt_type,complexity
    if(a):
        m=type_choice.find('Y',0)
    else:
        m=type_choice.find('X',0)
    if(m>=0):
        q=numpy.random.choice(grp_type,None,True,wt_type)  
        if(a==1 and q=="[H]"):
            q=grp_type[1]
        type_choice=type_choice[0:m]+q+type_choice[m+1:]  
        if(q!="[H]"):
            complexity+=1
        return 1
    return 0

         
    
def select_core_grp(s) :
    global core_grp,core_type
    n=len(core_grp)
    indx = [i for i in range(n) if core_grp[i]==s]
    i=random.choice(indx)
    q=core_type[i]
    return q

def select_complexity(x):
    global complexity_max,beta
    x0=3
    if(x<=x0):
        return 1
    dx=(x-x0)/(complexity_max-x0)
    if(x<=complexity_max):
        if(random.random()<math.exp(-beta*dx)):
            return 1
    return 0

def accept_mol(smiles,k):
    global  mol_set,complexity_set
    global natom_set,n_complexity
    global nring
    
    m = Chem.MolFromSmiles(smiles,sanitize=False)
    if m is None:
        return 0
    else:
        try:
            Chem.SanitizeMol(m)
            smiles=Chem.MolToSmiles(m)
            n = m.GetNumAtoms() 
            rings=m.GetRingInfo()
            for r in rings.AtomRings():
                if(len(r) not in nring):  # exclude some rings
                    return 0
            if smiles not in mol_set: 
                print("(",k+1,")",complexity,smiles)
                mol_set.append(smiles)
                complexity_set.append(complexity)
                natom_set.append(n)
                n_complexity[complexity]+=1
                return 1   
            else:
                return 0
        except:
             return 0

    
    
def recognize_key(content,key): 
    n=len(content)
    for i in range(n):
        columns=return_columns(content,i)
        if columns[0].count(key):
            return i
    return -1  

def return_columns(content,i):
         line=content[i] 
         columns=line.rstrip().split(',')
         return columns
     
def return_column(content,i,n):
         line=content[i] 
         columns=line.rstrip().split(',')
         return columns[n]
     
def label_cyclic_type(grp,z): 
    global nring    
    atom=["C","N","O"]
    bad="+-123456789"
    atoms=[]
    n=len(grp)
    m=grp.find("*",0)
    
    if(m<0):
         for i in range(n):
             x=grp[i]
             x.upper()    
             if x in atom:
                 if(i<n-1):
                     x1=grp[i+1]
                 else:
                     x1=" "
                 if x1 not in bad:
                    atoms.append(i)
    
         nat=len(atoms)
         if(nat<min(nring)-2):
             return ""
         i=random.randint(2,nat-1) 
         m=atoms[i]
         m1=m+1
    else:
        m1=m+1
        m-=1
        
    if m<n-1:
        s=grp[0:m+1]+z+grp[m1:]
    else:
        s=grp+z

    return s
    
    

def find_pos(t):
    global unique_grp,pos_seq,grp_seq,t_seq
    pos_seq=[]
    grp_seq=[]
    t_seq=[]
    for s in unique_grp:
        u=0
        while(1):
            u=t.find("("+s+")",u)
            if(u>=0):
                pos_seq.append(u)
                grp_seq.append(s)
                u+=1
            else:
                break
    pos_seq,grp_seq=zip(*sorted(zip(pos_seq,grp_seq)))
    
    m=0
    for i in range(len(pos_seq)):
        n=pos_seq[i]
        t_seq.append(t[m:n])
        m=n+2+len(grp_seq[i])
    if(m<=len(t)):
        t_seq.append(t[m:])
    
    
def one_cyclization(t,a,b):
    global layer_seq,n_seq
    s="%"+str(10+n_seq)
    if ("%" not in layer_seq[a]) and ("%" not in layer_seq[b]):
        layer_seq[a]="("+grp_seq[a]+s+")"
        layer_seq[b]=s
        n_seq+=1

def do_cyclization(t):
    global cyc_a,cyc_b,freq_cyc,cyc_index
    global n_seq,grp_seq,layer_seq
    global layer
    layer=t
    
    if(len(cyc_index)==0):
        return 0
    
    layer_seq=list(grp_seq)
    n_seq=0
    n=len(pos_seq)
    
    for i in range(n):
         layer_seq[i]="("+layer_seq[i]+")"
    
    for i in range(len(freq_cyc)):
        j=numpy.random.choice(cyc_index,None,True,freq_cyc)
        a=min(cyc_a[j],cyc_b[j])
        b=max(cyc_a[j],cyc_b[j])
        if(a>0):
           one_cyclization(t,a-1,b-1)
            
    if(n_seq):
        layer=""
        for i in range(n):
            layer+=t_seq[i]+layer_seq[i]
        if(len(t_seq)>n):
              layer+=t_seq[n]
    else:
        return 0
    
    return n_seq
    

def main(): 
    
    global core_grp,core_type,unique_grp,grp_symm,grp_cast
    global wt_type,grp_type
    global core,layer,complexity
    global complexity_max,beta
    global mol_set,complexity_set,natom_set,n_complexity
    global cyc_index,cyc_a,cyc_b,freq_cyc
    global nring
    
    
    inp=input("list name: ")
    f=open(inp+'.csv','r')
    content=f.readlines()
    f.close()
    
    complexity_max=8
    beta=4
    nmax=25000
    nring=[5,6,7]
    name_root="BE"
    
    cwd=os.getcwd()
    
    mol_set=[]
    complexity_set=[]
    natom_set=[]
    n_complexity=[]
    
    core_grp=[]
    core_type=[]
    grp_type=[]
    wt_type=[]
    unique_grp=[]
    grp_symm=[]
    grp_cast=[]
     
  
    n=len(content)
    i=recognize_key(content,"core") 
    core=return_column(content,i,1)
    print("\ncore = ",core)

    k=0
    for i in range(recognize_key(content,"group")+1,n):
        columns=return_columns(content,i)
        u=columns[0]
        if(len(u)):
            core_grp.append(u)
            core_type.append(columns[1])
            if u not in unique_grp:
                unique_grp.append(u)
        else:
            break

    find_pos(core)
    
    for grp in unique_grp:
        print("\ngroup ",grp)
        print(core_grp)
        
    
    k=0
    wt=0
    for i in range(recognize_key(content,"types")+1,n):
         columns=return_columns(content,i)
         if(len(columns[0])):
             grp_type.append(columns[0])
             u=float(columns[1])
             wt_type.append(u)
             wt+=u
             k+=1 
         else:
            break
 
    print("\ntype\t\tweight")
    for i in range(k):
        wt_type[i]/=wt
        print("%s\t\t%.3f" % (grp_type[i],wt_type[i]))
    
    k=0
    print("\nsymmetries:\nposition\tgroup\trelation")
    for i in range(recognize_key(content,"position")+1,n): 
        columns=return_columns(content,i)
        if(len(columns[0])):
            grp_symm.append(int(columns[1]))
            grp_cast.append(columns[2])
            print("%d\t\t%s\t%d" % (k+1,grp_cast[k],grp_symm[k]))
            k+=1
        else:
            break
    
        
    k=0
    wt=0
    cyc_a=[]
    cyc_b=[]
    cyc_index=[]
    freq_cyc=[]
    for i in range(recognize_key(content,"cycle")+1,n): 
        columns=return_columns(content,i)
        if(len(columns[0])):
            cyc_a.append(int(columns[0]))
            cyc_b.append(int(columns[1]))
            u=float(columns[2])
            freq_cyc.append(u)
            wt+=u
            cyc_index.append(k)
            k+=1
        else:
            break

    if(k):
        print("\ncyclization:\nposition a\tposition b\tweight")
        for i in range(k):
            freq_cyc[i]/=wt
            print("%d\t\t%d\t\t%.3f" % (cyc_a[i],cyc_b[i],freq_cyc[i]))
    
    print("\n")
    
    k=0
    for i in range(complexity_max+1):
        n_complexity.append(0)
        
    while k<nmax:
        if(once()):
            if(select_complexity(complexity)):
              mp=pybel.readstring("smi",layer)
              smi=mp.write()
              k+=accept_mol(smi,k)
                    
    natom_set,complexity_set,mol_set=zip(*sorted(zip(natom_set,complexity_set,mol_set)))
            
    f=open('complexity_histogram.csv','w')
    print("score,frequency",file=f)
    for i in range(complexity_max+1):
        print("%d,%3e" % (i,float(n_complexity[i])/nmax), file=f)
    f.close()
    
    f=open('results.csv','w')
    print("name,smiles,complexity,atoms",file=f)
    k=len(mol_set)
    for i in range(k):
        name=name_root+str(i+1)
        smiles=mol_set[i]
        complexity=complexity_set[i]
        n = natom_set[i]
        print("%s,%s,%d,%d" %(name,smiles,complexity,n),file=f)
    f.close()

    
main()