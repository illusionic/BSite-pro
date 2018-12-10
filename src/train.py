# -*- coding: utf-8 -*-
"""
Created on Sun Nov 11 00:46:40 2018

@author: Illusion
"""
import re
from Bio import SeqIO
import csv
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
import numpy as np
import math,copy
from sklearn.metrics import make_scorer, accuracy_score, precision_score, recall_score, f1_score
import pickle

AALetter = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]


##################################################CTD#######################################################################
AALetter=["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]

_Hydrophobicity={'1':'RKEDQN','2':'GASTPHY','3':'CLVIMFW'}  #'1'stand for Polar; '2'stand for Neutral, '3' stand for Hydrophobicity

_NormalizedVDWV={'1':'GASTPD','2':'NVEQIL','3':'MHKFRYW'} #'1'stand for (0-2.78); '2'stand for (2.95-4.0), '3' stand for (4.03-8.08)

_Polarity={'1':'LIFWCMVY','2':'CPNVEQIL','3':'KMHFRYW'} #'1'stand for (4.9-6.2); '2'stand for (8.0-9.2), '3' stand for (10.4-13.0)

_Charge={'1':'KR','2':'ANCQGHILMFPSTWYV','3':'DE'} #'1'stand for Positive; '2'stand for Neutral, '3' stand for Negative

_SecondaryStr={'1':'EALMQKRH','2':'VIYCWFT','3':'GNPSD'} #'1'stand for Helix; '2'stand for Strand, '3' stand for coil

_SolventAccessibility={'1':'ALFCGIVW','2':'RKQEND','3':'MPSTHY'} #'1'stand for Buried; '2'stand for Exposed, '3' stand for Intermediate

_Polarizability={'1':'GASDT','2':'CPNVEQIL','3':'KMHFRYW'} #'1'stand for (0-0.108); '2'stand for (0.128-0.186), '3' stand for (0.219-0.409)


##You can continuely add other properties of AADs to compute descriptors of protein sequence.

_AATProperty=(_Hydrophobicity,_NormalizedVDWV,_Polarity,_Charge,_SecondaryStr,_SolventAccessibility,_Polarizability)

_AATPropertyName=('_Hydrophobicity','_NormalizedVDWV','_Polarity','_Charge','_SecondaryStr','_SolventAccessibility','_Polarizability')




def StringtoNum(ProteinSequence,AAProperty): 
	"""
	Tranform the protein sequence into the string form such as 32123223132121123.
	AAProperty is a dict form containing classifciation of amino acids such as _Polarizability. result is a string such as 123321222132111123222
	"""
	
	hardProteinSequence=copy.deepcopy(ProteinSequence)
	for k,m in AAProperty.items():
		for index in str(m):
			hardProteinSequence=str.replace(hardProteinSequence,index,k)
	TProteinSequence=hardProteinSequence

	return TProteinSequence

###################################################################################################################
def CalculateComposition(ProteinSequence,AAProperty,AAPName):
	"""
	A method used for computing composition descriptors.
	AAProperty is a dict form containing classifciation of amino acids such as _Polarizability.	
	AAPName is a string used for indicating a AAP name.	
    result is a dict form containing composition descriptors based on the given property.
	"""
	TProteinSequence=StringtoNum(ProteinSequence,AAProperty)
	Result={}
	Num=len(TProteinSequence)
	Result[AAPName+'C'+'1']=round(float(TProteinSequence.count('1'))/Num,3)
	Result[AAPName+'C'+'2']=round(float(TProteinSequence.count('2'))/Num,3)
	Result[AAPName+'C'+'3']=round(float(TProteinSequence.count('3'))/Num,3)
	return Result

def CalculateTransition(ProteinSequence,AAProperty,AAPName):  
	"""
	A method used for computing transition descriptors	
	AAProperty is a dict form containing classifciation of amino acids such as _Polarizability.
    
	AAPName is a string used for indicating a AAP name.
    
	Output:result is a dict form containing transition descriptors based on the given property.
	"""
	
	TProteinSequence=StringtoNum(ProteinSequence,AAProperty)
	Result={}
	Num=len(TProteinSequence)
	CTD=TProteinSequence
	Result[AAPName+'T'+'12']=round(float(CTD.count('12')+CTD.count('21'))/(Num-1),3)
	Result[AAPName+'T'+'13']=round(float(CTD.count('13')+CTD.count('31'))/(Num-1),3)
	Result[AAPName+'T'+'23']=round(float(CTD.count('23')+CTD.count('32'))/(Num-1),3)
	return Result



def CalculateDistribution(ProteinSequence,AAProperty,AAPName):  
	
	"""
	A method used for computing distribution descriptors.
	"""
	TProteinSequence=StringtoNum(ProteinSequence,AAProperty)
	Result={}
	Num=len(TProteinSequence)
	temp=('1','2','3')
	for i in temp:
		num=TProteinSequence.count(i)
		ink=1
		indexk=0
		cds=[]
		while ink<=num:
			indexk=str.find(TProteinSequence,i,indexk)+1
			cds.append(indexk)
			ink=ink+1
				
		if cds==[]:
			Result[AAPName+'D'+i+'001']=0
			Result[AAPName+'D'+i+'025']=0
			Result[AAPName+'D'+i+'050']=0
			Result[AAPName+'D'+i+'075']=0
			Result[AAPName+'D'+i+'100']=0
		else:
				
			Result[AAPName+'D'+i+'001']=round(float(cds[0])/Num*100,3)
			Result[AAPName+'D'+i+'025']=round(float(cds[int(math.floor(num*0.25))-1])/Num*100,3)
			Result[AAPName+'D'+i+'050']=round(float(cds[int(math.floor(num*0.5))-1])/Num*100,3)
			Result[AAPName+'D'+i+'075']=round(float(cds[int(math.floor(num*0.75))-1])/Num*100,3)
			Result[AAPName+'D'+i+'100']=round(float(cds[-1])/Num*100,3)

	return Result

####################################################################################################
def CalculateCompositionHydrophobicity(ProteinSequence):
	
	"""
	A method used for calculating composition descriptors based on Hydrophobicity of AADs.
	"""
	
	result=CalculateComposition(ProteinSequence,_Hydrophobicity,'_Hydrophobicity')
	return result
	
def CalculateCompositionNormalizedVDWV(ProteinSequence):
	"""
	A method used for calculating composition descriptors based on NormalizedVDWV of AADs.
	"""
	result=CalculateComposition(ProteinSequence,_NormalizedVDWV,'_NormalizedVDWV')
	return result
	
def CalculateCompositionPolarity(ProteinSequence):
	"""
	A method used for calculating composition descriptors based on Polarity of AADs.	
	"""
	
	result=CalculateComposition(ProteinSequence,_Polarity,'_Polarity')
	return result
	
def CalculateCompositionCharge(ProteinSequence):
	"""
	A method used for calculating composition descriptors based on Charge of AADs.

	"""
	
	result=CalculateComposition(ProteinSequence,_Charge,'_Charge')
	return result
	
def CalculateCompositionSecondaryStr(ProteinSequence):
	"""
	###############################################################################################
	A method used for calculating composition descriptors based on SecondaryStr of AADs.
	"""
	
	result=CalculateComposition(ProteinSequence,_SecondaryStr,'_SecondaryStr')
	return result
	
def CalculateCompositionSolventAccessibility(ProteinSequence):
	"""
	A method used for calculating composition descriptors based on SolventAccessibility of  AADs.	
	Output:result is a dict form containing Composition descriptors based on SolventAccessibility.
	"""
	
	result=CalculateComposition(ProteinSequence,_SolventAccessibility,'_SolventAccessibility')
	return result
##################################################################################################
def CalculateCompositionPolarizability(ProteinSequence):
	"""
	A method used for calculating composition descriptors based on Polarizability of AADs.
	"""
	
	result=CalculateComposition(ProteinSequence,_Polarizability,'_Polarizability')
	return result


#######################################################################################################
def CalculateTransitionHydrophobicity(ProteinSequence):
	"""
	A method used for calculating Transition descriptors based on Hydrophobicity ofAADs.
	"""
	
	result=CalculateTransition(ProteinSequence,_Hydrophobicity,'_Hydrophobicity')
	return result
	
def CalculateTransitionNormalizedVDWV(ProteinSequence):
	"""
	###############################################################################################
	A method used for calculating Transition descriptors based on NormalizedVDWV of AADs.
	"""
	
	result=CalculateTransition(ProteinSequence,_NormalizedVDWV,'_NormalizedVDWV')
	return result
	
def CalculateTransitionPolarity(ProteinSequence):
	"""
	A method used for calculating Transition descriptors based on Polarity of AADs.
	"""
	
	result=CalculateTransition(ProteinSequence,_Polarity,'_Polarity')
	return result
	
def CalculateTransitionCharge(ProteinSequence):
	"""
	A method used for calculating Transition descriptors based on Charge of AADs.
	"""
	
	result=CalculateTransition(ProteinSequence,_Charge,'_Charge')
	return result
	
def CalculateTransitionSecondaryStr(ProteinSequence):
	"""
	A method used for calculating Transition descriptors based on SecondaryStr of AADs.
	"""
	
	result=CalculateTransition(ProteinSequence,_SecondaryStr,'_SecondaryStr')
	return result
	
def CalculateTransitionSolventAccessibility(ProteinSequence):
	"""
	A method used for calculating Transition descriptors based on SolventAccessibility of  AADs.
	"""
	
	result=CalculateTransition(ProteinSequence,_SolventAccessibility,'_SolventAccessibility')
	return result
	
def CalculateTransitionPolarizability(ProteinSequence):
	"""
	A method used for calculating Transition descriptors based on Polarizability of AADs.
	"""
	
	result=CalculateTransition(ProteinSequence,_Polarizability,'_Polarizability')
	return result


############################################################################################################
def CalculateDistributionHydrophobicity(ProteinSequence):
	"""
	A method used for calculating Distribution descriptors based on Hydrophobicity of AADs.
	"""
	
	result=CalculateDistribution(ProteinSequence,_Hydrophobicity,'_Hydrophobicity')
	return result
	
def CalculateDistributionNormalizedVDWV(ProteinSequence):
	"""
	A method used for calculating Distribution descriptors based on NormalizedVDWV of AADs.
	"""
	
	result=CalculateDistribution(ProteinSequence,_NormalizedVDWV,'_NormalizedVDWV')
	return result
	
def CalculateDistributionPolarity(ProteinSequence):
	"""
	A method used for calculating Distribution descriptors based on Polarity of AADs.
	"""
	
	result=CalculateDistribution(ProteinSequence,_Polarity,'_Polarity')
	return result
	
def CalculateDistributionCharge(ProteinSequence):
	"""
	A method used for calculating Distribution descriptors based on Charge of ADDs.
	"""
	
	result=CalculateDistribution(ProteinSequence,_Charge,'_Charge')
	return result
	
def CalculateDistributionSecondaryStr(ProteinSequence):
	"""
	A method used for calculating Distribution descriptors based on SecondaryStr of AADs.
	"""
	
	result=CalculateDistribution(ProteinSequence,_SecondaryStr,'_SecondaryStr')
	return result
	
def CalculateDistributionSolventAccessibility(ProteinSequence):
	
	"""
	A method used for calculating Distribution descriptors based on SolventAccessibility of  AADs.
	"""
	
	result=CalculateDistribution(ProteinSequence,_SolventAccessibility,'_SolventAccessibility')
	return result
	
def CalculateDistributionPolarizability(ProteinSequence):
	"""
	A method used for calculating Distribution descriptors based on Polarizability of AADs.
	"""
	
	result=CalculateDistribution(ProteinSequence,_Polarizability,'_Polarizability')
	return result

##################################################################################################

def CalculateC(ProteinSequence):
	"""
	Calculate all composition descriptors based on seven different properties of AADs.
	"""
	result={}
	result.update(CalculateCompositionPolarizability(ProteinSequence))
	result.update(CalculateCompositionSolventAccessibility(ProteinSequence))
	result.update(CalculateCompositionSecondaryStr(ProteinSequence))
	result.update(CalculateCompositionCharge(ProteinSequence))
	result.update(CalculateCompositionPolarity(ProteinSequence))
	result.update(CalculateCompositionNormalizedVDWV(ProteinSequence))
	result.update(CalculateCompositionHydrophobicity(ProteinSequence))
	return result
	
def CalculateT(ProteinSequence):
	"""
	Calculate all transition descriptors based on seven different properties of AADs.
	"""
	result={}
	result.update(CalculateTransitionPolarizability(ProteinSequence))
	result.update(CalculateTransitionSolventAccessibility(ProteinSequence))
	result.update(CalculateTransitionSecondaryStr(ProteinSequence))
	result.update(CalculateTransitionCharge(ProteinSequence))
	result.update(CalculateTransitionPolarity(ProteinSequence))
	result.update(CalculateTransitionNormalizedVDWV(ProteinSequence))
	result.update(CalculateTransitionHydrophobicity(ProteinSequence))
	return result
	
def CalculateD(ProteinSequence):
	"""
	Calculate all distribution descriptors based on seven different properties of AADs.
	"""
	result={}
	result.update(CalculateDistributionPolarizability(ProteinSequence))
	result.update(CalculateDistributionSolventAccessibility(ProteinSequence))
	result.update(CalculateDistributionSecondaryStr(ProteinSequence))
	result.update(CalculateDistributionCharge(ProteinSequence))
	result.update(CalculateDistributionPolarity(ProteinSequence))
	result.update(CalculateDistributionNormalizedVDWV(ProteinSequence))
	result.update(CalculateDistributionHydrophobicity(ProteinSequence))
	return result


def CalculateCTD(ProteinSequence):
	"""
	Calculate all CTD descriptors based on seven different properties of AADs.
	"""
	result={}
	result.update(CalculateCompositionPolarizability(ProteinSequence))
	result.update(CalculateCompositionSolventAccessibility(ProteinSequence))
	result.update(CalculateCompositionSecondaryStr(ProteinSequence))
	result.update(CalculateCompositionCharge(ProteinSequence))
	result.update(CalculateCompositionPolarity(ProteinSequence))
	result.update(CalculateCompositionNormalizedVDWV(ProteinSequence))
	result.update(CalculateCompositionHydrophobicity(ProteinSequence))
	result.update(CalculateTransitionPolarizability(ProteinSequence))
	result.update(CalculateTransitionSolventAccessibility(ProteinSequence))
	result.update(CalculateTransitionSecondaryStr(ProteinSequence))
	result.update(CalculateTransitionCharge(ProteinSequence))
	result.update(CalculateTransitionPolarity(ProteinSequence))
	result.update(CalculateTransitionNormalizedVDWV(ProteinSequence))
	result.update(CalculateTransitionHydrophobicity(ProteinSequence))
	result.update(CalculateDistributionPolarizability(ProteinSequence))
	result.update(CalculateDistributionSolventAccessibility(ProteinSequence))
	result.update(CalculateDistributionSecondaryStr(ProteinSequence))
	result.update(CalculateDistributionCharge(ProteinSequence))
	result.update(CalculateDistributionPolarity(ProteinSequence))
	result.update(CalculateDistributionNormalizedVDWV(ProteinSequence))
	result.update(CalculateDistributionHydrophobicity(ProteinSequence))
	return result



###############################################Conjoint Tried method #############################################

#a Dipole scale (Debye): -, Dipole<1.0; +, 1.0<Dipole<2.0; ++, 2.0<Dipole<3.0; +++, Dipole>3.0; +'+'+', Dipole>3.0 with opposite orientation.
#b Volume scale (Ã…3): -, Volume<50; +, Volume> 50.
#c Cys is separated from class 3 because of its ability to form disulfide bonds.
 
_repmat={1:["A",'G','V'],2:['I','L','F','P'],3:['Y','M','T','S'],4:['H','N','Q','W'],5:['R','K'],6:['D','E'],7:['C']}
def _Str2Num(proteinsequence):
	"""
	translate the amino acid letter into the corresponding class based on the
	
	given form.
	
	"""
	repmat={}
	for i in _repmat:
		for j in _repmat[i]:
			repmat[j]=i
			
	res=proteinsequence
	for i in repmat:
		res=res.replace(i,str(repmat[i]))
	return res


###############################################################################
def CalculateConjointTriad(proteinsequence):

    res={}
    proteinnum=_Str2Num(proteinsequence)
    for i in range(8):
        for j in range(8):
            for k in range(8):
                temp=str(i)+str(j)+str(k)
                res[temp]=proteinnum.count(temp)
    return res




###########################################################AAC############################################################################
def CalculateAAComposition(ProteinSequence):

    LengthSequence = len(ProteinSequence)
    Result = {}
    for i in AALetter:
        Result[i] = round((float(ProteinSequence.count(i)) / LengthSequence)*100, 2)
    return Result

#20 features


#############################################################################################
def CalculateDipeptideComposition(ProteinSequence):

    #400 dipepdite features

    LengthSequence = len(ProteinSequence)
    Result = {}
    for i in AALetter:
        for j in AALetter:
            Dipeptide = i + j
            Result[Dipeptide] = round((float(ProteinSequence.count(Dipeptide)) / (LengthSequence - 1))*100, 2)
    return Result



#############################################################################################

def Getkmers():

    #8000 tripeptide name
    kmers = list()
    for i in AALetter:
        for j in AALetter:
            for k in AALetter:
                kmers.append(i + j + k)
    return kmers



#############################################################################################
def GetSpectrumDict(proteinsequence):
  #8000 tripeptide features
     
    result = {}
    kmers = Getkmers()
    for i in kmers:
        result[i] = len(re.findall(i, proteinsequence)*100)
    return result



#############################################################################################
def CalculateAADipeptideComposition(ProteinSequence):

    #8420 all uni,bi,tri features

    result = {}
    result.update(CalculateAAComposition(ProteinSequence))
    result.update(CalculateDipeptideComposition(ProteinSequence))

    result.update(GetSpectrumDict(ProteinSequence))

    return result
  
#############################################################################################
#fbr=open('binding.txt','r')
#ff=open('ready_binding.txt','w')
    

####specify known labels file which contains protein id and interpro labels#################
fflbl=open('lbl.txt','r')
###creates a new file with same protein ids repeated with different associated IPR labels#################

def ready_file():
 
 for line in fflbl:
     
    c=''
    a=line
    splitted = a.split('\t')
    first = splitted[0]

    f=''
    fs = line[len(first)+1:(len(line)-2)]
    fs=list(fs.split(';'))

    for i in fs:
       if i!='':         
         with open('binding.txt','r') as fbr:
            for bind in fbr:
             if i in bind:
                aa=line
                splt = aa.split()
                f = splt[0]
                c=c+i+';'

                if f!='':
                 with open('ready_binding_all_dataset_updated.txt', 'a') as ff: 
                    ff.write(f+' '+i+';'+'\n')

######feature extraction given the input protein sequence e:g '013033' in .fasta format#######     
def extract_feature_for_prediction():

     
     for record in SeqIO.parse ( "predict_unknowns.fasta", "fasta" ):
         if 'O13033' in record.description:
            
            p=str(record.seq)
            CalculateAAComposition(p)
            CalculateDipeptideComposition(p)
            GetSpectrumDict(p)
            res = CalculateAADipeptideComposition(p )


            Ctried=CalculateConjointTriad(protein)
            CTD=CalculateCTD(protein)

            
            
          
            flist=list(res.values())+list(Ctried.values())+list(CTD.values())
           
            
        
            a=np.array(flist)
          

     return (a.reshape(1,-1))
                 
                 
                 
fflbl=open('ready_binding_all_dataset.txt','r') 

def store_feature_value():

 for line in fflbl:
     
     a=line
     splitted = a.split()
     first = splitted[0]
     
     for record in SeqIO.parse ( "example.fasta", "fasta" ):
         if line[:len(first)] in record.description:
            
            p=str(record.seq)
            CalculateAAComposition(p)
            CalculateDipeptideComposition(p)
            GetSpectrumDict(p)
            res = CalculateAADipeptideComposition(p )

            Ctried=CalculateConjointTriad(protein)
            CTD=CalculateCTD(protein)

            flist=list(res.values())+list(Ctried.values())+list(CTD.values())
           
            

            #fs = line[len(first)+1:(len(line)-2)].split(';')
            
            fs = line[len(first)+1:(len(line)-2)]
            flist.append(fs)
            flist.append(line[:len(first)])
            
            myFile = open('binding_feature_all_dataset_new.csv', 'a')
            with myFile:
                writer = csv.writer(myFile)
                writer.writerow(flist)

###################################################################################################
        
#store features name 
def store_feature_name():

 CalculateAAComposition(protein)
 CalculateDipeptideComposition(protein)

 a=[]    
 GetSpectrumDict(protein)

 AAC = CalculateAADipeptideComposition(protein)
 
 Ctried=CalculateConjointTriad(protein)
 
 CTD=CalculateCTD(protein)
 



 b=list(AAC.keys())

 d=list(Ctried.keys())
 e=list(CTD.keys()) 
 a=b+d+e

 
 a.append('labels')
 a.append('pid')
 
 
 myFile = open('binding_feature_all_dataset_new.csv', 'a')
 with myFile:
  writer = csv.writer(myFile)
  writer.writerow(a)

###################################################################################################
# end store feature name 
  
def predict_model():
 modelname = 'binding_site_trained_model.sav'   
 loaded_model = pickle.load(open(modelname, 'rb'))
 
 ynew = loaded_model.predict(extract_feature_for_prediction())
 print("your binding site against this 'O13033' protein id is ",ynew)

    

def train_model():
 import pandas as pd
 df = pd.read_csv('binding_feature_all_dataset_new.csv')
 #print(df)
 data = df.drop('pid', axis=1)
 
 data = data.drop('labels', axis=1)
 
 target_labels=df['labels']
 
 unique_target_labels=df['labels'].unique().tolist()
 
 print('unique labels',len(unique_target_labels))
 
 '''
 from sklearn.preprocessing import label_binarize
 target_labels = label_binarize(target_labels,target_labels)#[0, 1, 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42])
 '''
 #total_target_labels = target_labels.shape[1]

 data_train, data_test, target_labels_train, target_labels_test = train_test_split(data, target_labels,random_state=0,test_size=0.3)

 random_forest = RandomForestClassifier(n_jobs=2,n_estimators=600, random_state=0,max_features='auto',criterion='entropy',min_samples_leaf=2,oob_score=True)

 random_forest.fit(data_train,target_labels_train ) #y train_y.values.ravel())




 y_predict = random_forest.predict(data_test)


 print('precission: ',precision_score(target_labels_test, y_predict, average='micro'))

 print('Recall: ',recall_score(target_labels_test,y_predict,average='micro'))

 print('New predcition Accuracy (test): ',(accuracy_score(target_labels_test, y_predict)*100))

 training_predict = random_forest.predict(data_train)
 print('Training Accuracy is:',(accuracy_score(target_labels_train, training_predict)*100))

 print('F1 score: ',f1_score(target_labels_test, y_predict, average='micro'))

 modelname = 'binding_site_trained_model.sav'
 pickle.dump(random_forest, open(modelname, 'wb'))
 

  
  

if __name__ == "__main__":
    
 protein = 'MNTDQQPYQGQTDYTQGPGNGQSQEQDYDQYGQPLYPSQADGYYDPNVAAGTEADMYGQQ'
 print('processing: ') 
 
 #ready_file()
 
###to store feature names and values in a csv file###
 #store_feature_name()
 #store_feature_value()
 
 train_model()
 predict_model()


###########################################Plot ROC curve############################################
'''
from sklearn.metrics import roc_curve, auc
#from scipy import interp
#from itertools import cycle

fpr = dict()
tpr = dict()
roc_auc = dict()

total_target_labels = target_labels.shape[1]

for i in range(total_target_labels):
    fpr[i], tpr[i], _ = roc_curve(target_labels_test[:, i], y_predict[:, i])
    roc_auc[i] = auc(fpr[i], tpr[i])

fpr["micro"], tpr["micro"], _ = roc_curve(target_labels_test.ravel(), y_predict.ravel())
roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

plt.figure()

lw = 2
plt.plot(fpr["micro"], tpr["micro"],label='micro-average ROC curve (area = {0:0.2f})' ''.format(roc_auc["micro"]),color='deeppink', linestyle=':', linewidth=2)
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC')
plt.legend(loc="lower right")

plt.show()

plt.savefig("rocBinding.pdf")
'''


print("successfuly: ")
#fbr.close()  
#ff.close()   
fflbl.close()
#fflabel.close()

 
 

    
    
 

