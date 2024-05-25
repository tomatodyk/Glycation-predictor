from PySide2.QtGui import QIcon
import tkinter as tk
from tkinter import filedialog
import openpyxl
import onnxruntime as ort
from PySide2.QtCore import QFile
from PySide2.QtUiTools import QUiLoader
from PySide2.QtWidgets import QApplication, QMainWindow, QPushButton, QPlainTextEdit, QMessageBox, QComboBox
import math
import os
import numpy as np
import pandas as pd



def read_file(filepath):
    try:
        fp = open(filepath)
    except IOError:
        exit()
    else:
        fp=open(filepath)
        lines = fp.read().splitlines()
        sequence={}
        protein_name=''
        namelist=[]
        for line in lines:
            if line[0]=='>':
                if protein_name !='':
                    sequence[protein_name]=seq
                    namelist.append(protein_name)
                seq=''
                protein_name=line
            else:
                seq+=line
        sequence[protein_name]=seq
        namelist.append(protein_name)
        return sequence,namelist

def get_frag(protein_seq,size_windows,protein_name):
    frag_data=[]
    L=int((size_windows-1)/2)
    end=len(protein_seq)
    position=[AA for AA,v in enumerate(protein_seq) if v=='K']
    for i in range(len(position)):
        pos=position[i]
        if pos-L>=0 and end-pos-L>0:
            frag_data.append('>K'+str(pos+1))
            frag_data.append(protein_seq[pos-L:pos+L+1])
        if pos-L>=0 and end-pos-L<=0:
            sup_right='O'*(L+1-(end-pos))
            frag_data.append('>K'+str(pos+1))
            frag_data.append(protein_seq[pos-L:]+sup_right)
        if pos-L<0 and end-pos-L>0:
            sup_left='O'*(L-pos)
            frag_data.append('>K'+str(pos+1))
            frag_data.append(sup_left+protein_seq[:pos+L+1])
        if pos-L<0 and end-pos-L<=0:
            sup_left='O'*(L-pos)
            sup_right='O'*(L+1-(end-pos))
            frag_data.append('>K'+str(pos+1))
            frag_data.append(sup_left+protein_seq[:]+sup_right)
    return frag_data

def extract_predict(filepath,size_windows=31):
    sequence,namelist=read_file(filepath)
    frag=[]
    for protein_name in namelist:
        protein_seq=sequence[protein_name]
        frag_data=get_frag(protein_seq,size_windows,protein_name)
        frag+=frag_data
    return frag



# ---- AAC、PWAA、EBGW三种特征提取方式  ---- #
def AAC(frag):
    lines = frag
    L=len(lines[1])
    n=int(len(lines))
    AAs='ACDEFGHIKLMNPQRSTVWYO'
    m=len(AAs)
    aac=np.zeros((n,m))
    for i in range(n):
        for j in range(m):
            frequency=lines[i].count(AAs[j])
            frequency=float('%.2f'%frequency)
            aac[i][j]=frequency/L
    aac=aac[:,0:20]
    return np.array(aac)

def PWAA(frag):
    lines = frag
    L=len(lines[1])
    n=int(len(lines))
    AAs='ACDEFGHIKLMNPQRSTVWYO'
    l=int((L-1)/2)
    data=np.zeros((n,21))
    for i in range(n):
        for k in range(len(AAs)):
            pos=[ii for ii,v in enumerate(lines[i]) if v==AAs[k]]
            pos2=[jj+1 for jj in pos]
            p=[]
            c=[]
            for j in pos2:
                if j>=1 and j<=l:
                    p.append(j-l-1)
                if j>l and j<=L:
                    p.append(j-l-1)
            for m in p:
                if m>=-l and m<=l:
                    S1=float('%.2f'%abs(m))
                    c.append(m+S1/l)
            S2=float('%.2f'%sum(c))
            data[i][k]=S2/(l*(l+1))
    return data
def EBGW(frag):
    lines = frag
    L = len(lines[1])
    l = L
    n = int(len(lines))
    C1 = 'AFGILMPVW'
    C2 = 'CNQSTY'
    C3 = 'HKR'
    C4 = 'DE'
    ucidata = []
    EBGW = []
    ucidata1 = []
    ucidata2 = []
    ucidata3 = []
    for i in range(n):
        ucida = []
        for j in range(l):
            pos = [ii for ii, v in enumerate(C1) if v == lines[i ][j]]
            pos1 = [ii for ii, v in enumerate(C2) if v == lines[i ][j]]
            pos2 = [ii for ii, v in enumerate(C3) if v == lines[i ][j]]
            pos3 = [ii for ii, v in enumerate(C4) if v == lines[i ][j]]
            if len(pos) == 1 or len(pos1) == 1:
                ucida.append(1)
            elif len(pos2) == 1 or len(pos3) == 1:
                ucida.append(0)
            else:
                ucida.append(0)
            pos = []
            pos1 = []
            pos2 = []
            pos3 = []

        ucidata1.append(ucida)

    for i in range(n):
        ucida1 = []
        for j in range(l):
            pos = [ii for ii, v in enumerate(C1) if v == lines[i][j]]
            pos1 = [ii for ii, v in enumerate(C3) if v == lines[ i ][j]]
            pos2 = [ii for ii, v in enumerate(C2) if v == lines[i ][j]]
            pos3 = [ii for ii, v in enumerate(C4) if v == lines[i][j]]
            if len(pos) == 1 or len(pos1) == 1:
                ucida1.append(1)
            elif len(pos2) == 1 or len(pos3) == 1:
                ucida1.append(0)
            else:
                ucida1.append(0)
            pos = []
            pos1 = []
            pos2 = []
            pos3 = []
        ucidata2.append(ucida1)
    for i in range(n):
        ucida2 = []
        for j in range(l):
            pos = [ii for ii, v in enumerate(C1) if v == lines[i][j]]
            pos1 = [ii for ii, v in enumerate(C3) if v == lines[ i ][j]]
            pos2 = [ii for ii, v in enumerate(C2) if v == lines[i ][j]]
            pos3 = [ii for ii, v in enumerate(C4) if v == lines[i][j]]
            if len(pos) == 1 or len(pos1) == 1:
                ucida2.append(1)
            elif len(pos2) == 1 or len(pos3) == 1:
                ucida2.append(0)
            else:
                ucida2.append(0)
            pos = []
            pos1 = []
            pos2 = []
            pos3 = []
        ucidata3.append(ucida2)
    ucidata = np.hstack((ucidata1, ucidata2, ucidata3))

    ur, uc = np.shape(ucidata)
    k1 = 5
    x1 = []
    x2 = []
    x3 = []
    for i in range(ur):
        x11 = []
        x22 = []
        x33 = []
        a = 0
        b = 0
        c = 0
        for j in range(int(k1)):
            a = sum(ucidata1[i][0:int(math.floor(l * (j + 1) / k1))]) / math.floor(l * (j + 1) / k1)
            b = sum(ucidata2[i][0:int(math.floor(l * (j + 1) / k1))]) / math.floor(l * (j + 1) / k1)
            c = sum(ucidata3[i][0:int(math.floor(l * (j + 1) / k1))]) / math.floor(l * (j + 1) / k1)
            x11.append(a)
            x22.append(b)
            x33.append(c)
        x1.append(x11)
        x2.append(x22)
        x3.append(x33)
    EBGW = np.hstack((x1, x2, x3))
    EBGW = np.array(EBGW)
    return EBGW

def read_file1(filepath):
        fp=open(filepath)
        lines = fp.read().splitlines()
        sequence=[]
        for line in lines:
            if line[0]!='>':
                sequence.append(line)
        return sequence

# ----特征选择 ---- #
def select(feature):
    pospath = 'index.csv'
    df = pd.read_csv(pospath)
    pos = df.values
    posL = []
    for i in range(len(pos)):
        posL.append(pos[i][0])
    new_feature = feature[:, posL]
    return new_feature

# ----特征提取 ---- #
def extract(inputfile):
    os.system('python ./Pse-in-One-2.0/nac.py ./'+inputfile+' Protein Kmer -k 2 -f tab -labels 0 -out PCBkmer.txt')
    df = pd.DataFrame(pd.read_csv('PCBkmer.txt', delimiter='\t',header=None))
    Kmer=df.values
    print (Kmer.shape)
    os.system('python ./Pse-in-One-2.0/nac.py ./'+inputfile+' Protein DR -max_dis 3 -f tab -labels 0 -out PCBDR.txt')
    df = pd.DataFrame(pd.read_csv('PCBDR.txt', delimiter='\t', header=None))
    DR = df.values
    print (DR.shape)
    frag = read_file1(inputfile)
    acc=AAC(frag)
    print(acc.shape)
    pwaa=PWAA(frag)
    print(pwaa.shape)
    ebgw=EBGW(frag)
    print(ebgw.shape)
    features=np.concatenate((acc,pwaa,ebgw,Kmer,DR),axis=1)
    print(features.shape)
    features=select(features)
    X = pd.DataFrame(features)
    X.to_csv("middle.csv" ,index=False, index_label=None, header=None)

# ---State类---- #
class State:
    namelist = []  # ----序号----#
    predlist = []  # ---保存预测概率
    y = []  # ---保存预测结果
    N=0
    def __init__(self):
        qfile_states = QFile("predict.ui")
        qfile_states.open(QFile.ReadOnly)
        qfile_states.close()
        self.ui = QUiLoader().load(qfile_states)
        self.ui.setObjectName("MainWindow")
        self.ui.setStyleSheet("#MainWindow{border-image:url(dna.jpg)}")#设置背景图片
        self.ui.button2.clicked.connect(self.load)
        self.ui.button3.clicked.connect(self.export)

    def load(self):
        # ----加载训练好的深度学习模型----#
        ort_session = ort.InferenceSession("model_TEST.onnx")
        # ----获取指定的文件名称----#
        root = tk.Tk()
        root.withdraw()
        # 打开文件对话框
        file_path = filedialog.askopenfilename(parent=root, title='选择文件')
        # 检查文件路径是否为空
        if file_path=="":
            print('未选择文件')
        else:
            frag = extract_predict(file_path, 31)
            newFile0 = open('sequence.txt', 'w')
            for i in range(0, len(frag), 2):
                self.namelist.append(frag[i])
            for i in range(0, len(frag)):
                newFile0.write("".join(frag[i]) + '\n')
            newFile0.close()
            # -----特征提取------#
            IN = 'sequence.txt'
            extract(IN)
            # ------读取特征------#
            df = pd.DataFrame(pd.read_csv('middle.csv', header=None))
            X_test = df.values
            # ------进行预测------  #
            self.N = len(X_test)
            a = []
            self.ui.comtext.clear()
            self.ui.comtext.appendPlainText('序号'.ljust(10) + '预测概率'.ljust(10) + '预测结果')
            for i in range(0, len(X_test)):
                input_data = X_test[i].astype(np.float32)
                P_test = [[input_data]]
                output = ort_session.run(None, {"INPUT": P_test})
                if output[0][0][0] > 0.5:
                    self.y.append('糖化')
                    self.predlist.append(output[0][0][0])
                    a.append(str(self.namelist[i]).ljust(10) + str(output[0][0][0]).ljust(17) + '糖化')
                else:
                    self.y.append('非糖化')
                    self.predlist.append(output[0][0][0])
                    a.append(str(self.namelist[i]).ljust(10) + str(output[0][0][0]).ljust(17) + '非糖化')
                self.ui.comtext.appendPlainText(str(a[i]))

    def export(self):
        files = [('*xlsx', '*.xlsx')]
        file = filedialog.asksaveasfile(filetypes=files,defaultextension=files)
        print(file.name)
        workbook=openpyxl.Workbook()
        sheet=workbook.active
        sheet['A1']='序号'
        sheet['B1']='预测概率'
        sheet['C1']='预测结果'
        for i in range(0, self.N):
            sheet.append([self.namelist[i],self.predlist[i],self.y[i]])
        workbook.save(file.name)
# ----主函数---- #
if __name__ == '__main__':
    app = QApplication([])
    app.setWindowIcon(QIcon('pro.jpg'))
    state = State()
    state.ui.show()
    app.exec_()
