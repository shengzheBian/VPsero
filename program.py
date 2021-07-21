#main program: recognize input, use function.
#!/usr/bin/env python
# coding=utf-8

#change log:2021-7-21: grep hldD not uniprot id to query border gene, because the difference of prokka version. 



#main program: recognize input, use function.


import argparse
import os
import glob
import time
import re
import pandas as pd
import numpy as np
import sys
import logging
###############ingore warning info
pd.set_option('mode.chained_assignment',None)
#parameter input and recognize
parser = argparse.ArgumentParser(description="Prediction for Serotype of Vibrio parahaemolyticus")

parser.add_argument("-i","--in_dir", metavar="<dir>",
                                 required=False,action="store",dest="in_dir",
                                 help="a director that contains all genome assemble fasta file,the work dir must be in software")

parser.add_argument("-p","--prokka_result", metavar="<dir>",
                                                  required=False,action="store",dest="prokka_dir",                                                              help="set a director that contains all prokka results, and program will skip the prokka step, the work dir must be in software")

parser.add_argument("-o","--out_dir", metavar="<dir>",
                                  required=True,action="store",dest="out_dir",
                                 help="a direcotr that generate analyze result")

parser.add_argument("-n","--thread_num", metavar="<int>",
                                  required=False,action="store",dest="thread_num",                                                              help="set the thread number when genome annotate, default is 1")


args = parser.parse_args()
in_dir = args.in_dir
prokka_dir = args.prokka_dir
o_dir= args.out_dir
thread_num = args.thread_num


#####################################normalize user input
cwd = os.getcwd()

if in_dir != None:
    if in_dir[0] !="/":# -i transfer to absolute dir
     in_dir = cwd+"/"+in_dir

    if in_dir[-1] != "/":
     in_dir = in_dir+"/"

if o_dir[0] != "/": # -o transfer to absolute dir
    o_dir = cwd+"/" + o_dir

if o_dir[-1] != "/":
    o_dir = o_dir+"/"

if prokka_dir != None:
    if prokka_dir[0] !="/":# -i transfer to absolute dir
        prokka_dir = cwd+"/"+prokka_dir
    
    if prokka_dir[-1] != "/":
        prokka_dir = prokka_dir+"/"
 

os.system("mkdir -p "+o_dir+"serotype_predict")

o_dir_all = o_dir+"serotype_predict/"


if thread_num == None:
    thread_num = 1
else:
    thread_num = int(thread_num)

####set logging
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
   #create a handler,write in log
log_name=o_dir_all+"VPsero.log"
fh = logging.FileHandler(log_name,mode="w")
fh.setLevel(logging.DEBUG)
   #define the log output format
#fmt = "%(asctime)-15s %(levelname)s %(filename)s %(lineno)d %(message)s"
fmt = "%(asctime)-15s %(levelname)s %(lineno)d %(message)s"
#datefmt = "%a %d %b %Y %H:%M:%S"
datefmt = "%H:%M:%S"
formatter = logging.Formatter(fmt,datefmt)
   #add logger in handler
fh.setFormatter(formatter)
logger.addHandler(fh)

####program_dir and scripts source dir
program_py_dir=sys.argv[0]   #2021-7-16: auto find the program.py dir,and change the work dir
if program_py_dir == "program.py":
    program_dir = ""
if program_py_dir[0] == "/":
    dir_split=program_py_dir.split("/")
    program_dir = "/".join(dir_split[0:len(dir_split)-1])+"/"
logger.info("the program_dir is "+program_dir)


#input fasta name can't contain "ffn"
#xxxx
#xxx

#########################################################
#function: input -i , -o, and -n, generate a dir which contain pokka result of every genome.
#middle progress:
def genome_annote(in_dir=in_dir,o_dir_all=o_dir_all,thread_num=1):
    ### make director
    annote_dir = o_dir_all+"01.annote/"
    os.system("mkdir -p "+annote_dir)
    

    ### test if prokka is in environment 
    status = os.system("prokka -v")
    try:
        if (status != 0):
            raise RuntimeError
    except RuntimeError as e:
        logger.error("Prokka analyse is wrong ,Please check your software,for example, did you add it in environment?")
        os._exit(0)

    ### make prokka shell in one file
    fna_list = os.listdir(in_dir)
    
    os.system("mkdir -p "+o_dir_all+"middle_file/")

    for i in fna_list:
        prokka_cmd = "echo \"prokka --force --outdir %s --fast --quiet --prefix %s --gcode 11 --locustag %s %s\">>%s"%(annote_dir+i,i,i,in_dir+i,o_dir_all+"middle_file/prokka.sh")
        status = os.system(prokka_cmd)
        try:
            if (status != 0):
                raise RuntimeError
        except RuntimeError as e:
            logger.error("Prokka scripts generate wrong, please check it")
            os._exit(0)

    ### run prokka by thread_number setted.
    if len(fna_list)%thread_num == 0:
        l_set = len(fna_list)/thread_num
    else:
        l_set = len(fna_list)//thread_num+1
    
    python_dir = os.getcwd()#save the raw python work dir.
    os.chdir(o_dir_all+"middle_file/")#change the work dir to output dir.
    split_cmd = "split -l %s -a 4 -d %s prokka_thread" % (str(l_set),o_dir_all+"middle_file/prokka.sh")
    os.system(split_cmd)#linux run cmd by order.
    os.system("rm prokka.sh")

    logger.debug ("run cmd:"+split_cmd)

    prokka_thread_cmd_file = os.listdir(o_dir_all+"middle_file/")
    for i in prokka_thread_cmd_file:
        status = os.system("sh %s &"%(i))
        try:
            if (status != 0):
                raise RuntimeError
        except RuntimeError as e:
            logger.error("Prokka analyse is wrong ,Please check your software,for example, did you add it in environment?")
            os._exit(0)

        print "sh %s &"%(i)
    while True:
        aaaaaaaaaaaaaaa=1
        if len(glob.glob(annote_dir+"*/*ffn")) == len(fna_list):
            time.sleep(5)
            print glob.glob(annote_dir+"*/*ffn")
            print len(glob.glob(annote_dir+"*/*ffn"))
            break
    ##remove middle file
    os.system("rm ../middle_file/*")
    logger.info ("prokka module is OK!!!")
            
    os.chdir (python_dir)
    
############################################################################################
#function 1.1:
    #introduce: copy the prokka result prepared by user to 01.annote
def copy_prokka_result(prokka_dir=prokka_dir,o_dir_all=o_dir_all,thread_num=1):
    ### make director
    annote_dir = o_dir_all+"01.annote/"
    os.system("mkdir -p "+annote_dir)
    os.system("mkdir -p "+o_dir_all+"middle_file/")
    ### cp prokka result to 01.aanote/
    os.system("cp -r "+prokka_dir+"* "+annote_dir)
#############################################################################################
#function2:
    #introduce:input the prokka result dir, specgene fasta file,  and return the blast m8 file.
    #middle progress:1.make blast db
                   # 2.blast
                   # 3.m8 file coverage calculate,filter.first line add.
                   # 4.data toushi sheet.
def query_spec_gene(in_dir=o_dir_all+"01.annote/",in_spec=program_dir+"source/spec_gene_final.fasta"):
    #0.check the ffn if is 0:break

    #1.cat all genenome fasta to 1.
    cat_cmd = "cat %s*/*ffn>%s" % (in_dir , o_dir_all+"middle_file/all_genome.fasta")
    logger.info("")
    #logger.info("The command is as following:")
    logger.debug(cat_cmd)
    cat_stat = os.system (cat_cmd)#OK!! have test!
    try:
        if (cat_stat != 0):
            raise RuntimeError
    except RuntimeError as e:
        logger.error("cat * ffn file from prokka is error! please check your input data in 01.annot/, makfesure if they are generate from prokka or prepared, and one strain is one fold with .ffn postfix, the format please refer VPsero's exampledata/prokka_result/ ") 
        os._exit(0)     
    
    os.system ("mkdir -p %s02.blastn"%(o_dir_all))
    logger.debug(cat_cmd)
    
    #2.blast,calculate,filter,data toushi sheet
    
    sh_cmd = "python %s -idb %s -iq %s -o %s -iL %s" % (program_dir+"scripts_of/s12_mkdb_blastn.py",o_dir_all+"middle_file/all_genome.fasta", in_spec, o_dir_all+"02.blastn/blastn.m8", in_spec)
    logger.info("")
    #logger.info("The command is as following:")
    logger.debug(sh_cmd)
    if program_dir == "":
        status=os.system (sh_cmd)
    else:
        os.chdir(program_dir)  #2021-7-16
        status=os.system (sh_cmd)
        os.chdir(o_dir_all)    #2021-7-16
    try:
        if (status != 0):
            raise RuntimeError
    except RuntimeError as e:
        logger.error("blastn is Error! please Check it")
        os._exit(0)

    logger.info("blastn OK!")
    
###################################function3###################
   #function3:
   #introduction: input the prokka gff dir, output two table: O cluster info; K cluster info.
   #middle progress: 1.extract 3 border gene into 3 gff file.
                    #2.transfer 3 gff file to tsv table.
                    #3.merge two boder gene info table.(O and K)
def border_analyze(in_dir=o_dir_all+"01.annote/"):
    #cat all gff to one.
    os.system("mkdir -p %s" %(o_dir_all+"03.border_analyze"))#ouput dir of this step.
    
    def one_border_process(gene, in_dir=o_dir_all+"01.annote/"):
        #1.filter this border gene to one gff.
        grep_gff_cmd = "cat %s|grep ID |grep %s >%s" % (in_dir+"*/*gff",gene,o_dir_all+"03.border_analyze/"+gene)
        logger.info("")
        #print("The command is as following:")
        logger.debug(grep_gff_cmd)
        os.system(grep_gff_cmd)
        
        #2.trans gff to tsv.
        f = open (o_dir_all+"03.border_analyze/"+gene)
        prokka_list = [i.strip() for i in f]
        border_gene_list = []

        for line in prokka_list:
          line_split = line.split()
          contig=line_split[0]
          geneid = re.findall(r"locus_tag=(\S+);",line)
          uniprotid = re.findall(r"UniProtKB:(\w+);",line)
          pos1 = line_split[3]
          pos2 = line_split[4]
          direct = line_split[6]
          VPnumber = "_".join(geneid[0].split("_")[0:-1])
          
          if len(geneid)!=1 or len(uniprotid)!=1:
              print "regress extract no or more than 1"
          border_gene_list.append([VPnumber,contig,geneid[0],uniprotid[0],pos1,pos2,direct])    
         # print VPnumber,contig,geneid[0],uniprotid[0],pos1,pos2,direct
        f.close()
        
        border_gene_df = pd.DataFrame(border_gene_list)
      #  print border_gene_df
        return border_gene_df   
    #run function
    coaD = one_border_process(gene="coaD")#the border gene of O,coaD
    hldD = one_border_process(gene="hldD")#the border gene of O and K,hldD #2021-7-21 changed
    glpX = one_border_process(gene="P0A9C9")# the border gene of K,glpX
    #test = pd.DataFrame([["05-3133.fna","LRAI01000025.2","05-3133.fna_02345","B7V2S6","174654","175136","+"]])#######test
    #coaD = coaD.append(test,ignore_index=True)########test

    #merge to border table, and rm no on contig. #test! work well!
    coaD_hldD_df = pd.merge(coaD, hldD, left_on=[0],right_on=[0],how="outer")
    coaD_hldD_df = coaD_hldD_df[coaD_hldD_df["1_x"] == coaD_hldD_df["1_y"] ]
    #print coaD_hldD_df

    hldD_glpX_df = pd.merge(hldD, glpX, left_on=[0],right_on=[0],how="outer")
    hldD_glpX_df = hldD_glpX_df[hldD_glpX_df["1_x"] == hldD_glpX_df["1_y"]]
    #print hldD_glpX_df

    #left merge the whole strain number, and generate the final table.
    all_strain_list = os.listdir(in_dir)
    all_strain_df = pd.DataFrame(all_strain_list)

    all_strain_O_df = pd.merge(all_strain_df, coaD_hldD_df, left_on=[0],right_on=[0],how="left")
    #print all_strain_O_df
    all_strain_K_df = pd.merge(all_strain_df, hldD_glpX_df, left_on=[0],right_on=[0],how="left")
    #print all_strain_K_df

    all_strain_O_df.to_excel(o_dir_all+"03.border_analyze/all_strain_O_info.xlsx")
    all_strain_K_df.to_excel(o_dir_all+"03.border_analyze/all_strain_K_info.xlsx")

#function 4
    #introduciton: rm the blastn result whose gene_id doesn't in CPS or LPS cluter or didn't extract whole gene cluster.
    #middle process: input the filter blastn result in 02, and filter the row according to the 03 gene cluster result.
            #1.read in three excel to df
            #2.filter blastn result(in K cluster or O cluster)
            #3.long blastn to pivot
            #4.generate the final result.
def combine_02_03(in_blastn = o_dir_all+"02.blastn/filter_blastn.xlsx",in_O_info = o_dir_all+"03.border_analyze/all_strain_O_info.xlsx", in_K_info = o_dir_all+"03.border_analyze/all_strain_K_info.xlsx"):
    #1.read in three excel
    path_1=unicode(in_blastn,"utf8")
    filter_blastn_df = pd.read_excel(path_1, sheet_name = "Sheet1")#have test, OK!
    filter_blastn_df["subject_strain_id"] = filter_blastn_df["subject_id"].apply(lambda x: "_".join(x.split("_")[0:-1]))
    #print filter_blastn_df

    path_2 = unicode(in_O_info,"utf8")
    O_cluter_df = pd.read_excel(path_2, sheet_name = "Sheet1")#have test, OK!
    O_cluter_df.columns = ["O_strain_id","O_coaD_contig","O_coaD_geneid","O_coaD_uniprotid","O_coaD_pos1","O_coaD_pos2","O_coaD_direct","O_hldD_contig","O_hldD_geneid","O_hldD_uniprotid","O_hldD_pos1","O_hldD_pos2","O_hldD_direct"]
    #print O_cluter_df

    path_3 = unicode(in_K_info,"utf8")
    K_cluter_df = pd.read_excel(path_3, sheet_name = "Sheet1")#have test, OK!
    K_cluter_df.columns = ["K_strain_id","K_hldD_contig","K_hldD_geneid","K_hldD_uniprotid","K_hldD_pos1","K_hldD_pos2","K_hldD_direct","K_glpX_contig","K_glpX_geneid","K_glpX_uniprotid","K_glpX_pos1","K_glpX_pos2","K_glpX_direct"]
    #print K_cluter_df

    #2.filter blastn result(in K cluster or O cluster)
    filter_blastn_O_cluster = pd.merge(filter_blastn_df, O_cluter_df, left_on=["subject_strain_id"],right_on=["O_strain_id"],how="left")
    filter_blastn_OK_cluster = pd.merge(filter_blastn_O_cluster, K_cluter_df, left_on=["subject_strain_id"],right_on=["K_strain_id"],how="left")
    #print filter_blastn_OK_cluster
    #filter_blastn_OK_cluster.to_excel(o_dir_all+"03.border_analyze/filter_blastn_OK_cluster.xlsx")
    
    #2.A O predict result:
    blastn_O_df = filter_blastn_OK_cluster[filter_blastn_OK_cluster["query_id"].apply(lambda x:x[0])=="O"]
    blastn_O_dropna_df = blastn_O_df.dropna(axis=0,subset=["O_coaD_geneid"])#rm no O cluster OK!!!
    #print blastn_O_dropna_df
    blastn_O_dropna_df["subject_locus"] = blastn_O_dropna_df["subject_id"].apply(lambda x:int(x.split("_")[-1]))
    blastn_O_dropna_df["O_coaD_locus"] = blastn_O_dropna_df["O_coaD_geneid"].apply(lambda x:int(x.split("_")[-1]))
    blastn_O_dropna_df["O_hldD_locus"] = blastn_O_dropna_df["O_hldD_geneid"].apply(lambda x:int(x.split("_")[-1]))
    blastn_O_dropna_df["O_max"] = blastn_O_dropna_df[["O_coaD_locus","O_hldD_locus"]].max(axis = 1)
    blastn_O_dropna_df["O_min"] = blastn_O_dropna_df[["O_coaD_locus","O_hldD_locus"]].min(axis = 1)
    #print blastn_O_dropna_df
    blastn_O_dropna_df = blastn_O_dropna_df[(blastn_O_dropna_df["subject_locus"]<=blastn_O_dropna_df["O_max"]) & (blastn_O_dropna_df["subject_locus"]>=blastn_O_dropna_df["O_min"])]#rm the subject not in O cluster, have test, work well!!!
    #print blastn_O_dropna_df
    #blastn_O_dropna_df.to_excel(o_dir_all+"03.border_analyze/blastn_O_dropna_df.xlsx")

    #2.B K predict result:
    blastn_K_df = filter_blastn_OK_cluster[filter_blastn_OK_cluster["query_id"].apply(lambda x:x[0])=="K"]
    blastn_K_dropna_df = blastn_K_df.dropna(axis=0,subset=["K_glpX_geneid"])#rm no K cluster OK!!!
    #print blastn_K_dropna_df
    blastn_K_dropna_df["subject_locus"] = blastn_K_dropna_df["subject_id"].apply(lambda x:int(x.split("_")[-1]))
    blastn_K_dropna_df["K_hldD_locus"] = blastn_K_dropna_df["K_hldD_geneid"].apply(lambda x:int(x.split("_")[-1]))
    blastn_K_dropna_df["K_glpX_locus"] = blastn_K_dropna_df["K_glpX_geneid"].apply(lambda x:int(x.split("_")[-1]))
    blastn_K_dropna_df["K_max"] = blastn_K_dropna_df[["K_hldD_locus","K_glpX_locus"]].max(axis = 1)
    blastn_K_dropna_df["K_min"] = blastn_K_dropna_df[["K_hldD_locus","K_glpX_locus"]].min(axis = 1)
    #print blastn_K_dropna_df

    blastn_K_dropna_df = blastn_K_dropna_df[(blastn_K_dropna_df["subject_locus"]<=blastn_K_dropna_df["K_max"]) & (blastn_K_dropna_df["subject_locus"]>=blastn_K_dropna_df["K_min"])]#rm the subject not in K cluster ,have test, work well!!!!
    #print blastn_K_dropna_df

    #blastn_K_dropna_df.to_excel(o_dir_all+"03.border_analyze/blastn_K_dropna_df.xlsx")

     #3 long to pivot
       #3.A O
    blastn_O_pivot_df = blastn_O_dropna_df.pivot_table (index = "subject_strain_id", columns = "query_id", values = "subject_id", aggfunc = "count")    
    #print blastn_O_pivot_df
    print "----------------------------------"
    row_list = [a for a in blastn_O_pivot_df.index]
    row_num = len(blastn_O_pivot_df.index)
    final_predict_list = []
    for i in range(0,row_num):
        line = blastn_O_pivot_df.iloc[i]
        line_dropna = line.dropna()
        predict_list = [b.split("*")[0] for b in line_dropna.index]
        predict_result = ",".join(predict_list)
        out_strain_name = row_list[i]
        final_predict_list.append([out_strain_name,predict_result])
    #print final_predict_list
    O_predict_df = pd.DataFrame(final_predict_list,columns = ["strain_id","predict_O_result"])
    #print O_predict_df


       #3.B K
    blastn_K_pivot_df = blastn_K_dropna_df.pivot_table (index = "subject_strain_id", columns = "query_id", values = "subject_id", aggfunc = "count")
    #print blastn_K_pivot_df
    row_list2 = [d for d in blastn_K_pivot_df.index]
    row_num2 = len(row_list2)
    final_predict_list2 = []
    for e in range(0,row_num2):
        line2 = blastn_K_pivot_df.iloc[e]
        line2_dropna = line2.dropna()
        predict_list = [f.split("*")[0] for f in line2_dropna.index]
        predict_result = ",".join(predict_list)
        #print predict_result
        out_strain_name = row_list2[e]
        final_predict_list2.append([out_strain_name,predict_result])
    #print final_predict_list2
    K_predict_df = pd.DataFrame(final_predict_list2,columns = ["strain_id","predict_K_result"])
    #print K_predict_df

    #4.final output result
    all_strain_list = os.listdir(o_dir_all+"/01.annote")    
    #print all_strain_list
    all_stain_df = pd.DataFrame([[ele] for ele in all_strain_list],columns=["strain_name"])
    #print K_predict_df
    #print "-"*50
    merge_df1 = pd.merge(all_stain_df, O_cluter_df, left_on=["strain_name"],right_on=["O_strain_id"],how="left")
    merge_df2 = pd.merge(merge_df1,K_cluter_df,left_on=["strain_name"],right_on=["K_strain_id"],how="left")
    #print merge_df2
    
    merge_df3 =  pd.merge(merge_df2,O_predict_df,left_on=["strain_name"],right_on=["strain_id"],how="left")
    #print merge_df3

    merge_output = pd.merge(merge_df3,K_predict_df,left_on=["strain_name"],right_on=["strain_id"],how="left")
    #print merge_output
    #print merge_output.columns
    columns_list = [u'strain_name', u'O_coaD_contig', u'O_coaD_geneid',u'O_coaD_uniprotid', u'O_coaD_pos1', u'O_coaD_pos2', u'O_coaD_direct',
            u'O_hldD_contig', u'O_hldD_geneid', u'O_hldD_uniprotid', u'O_hldD_pos1',
            u'O_hldD_pos2', u'O_hldD_direct', u'K_hldD_contig',
            u'K_hldD_geneid', u'K_hldD_uniprotid', u'K_hldD_pos1', u'K_hldD_pos2',
            u'K_hldD_direct', u'K_glpX_contig', u'K_glpX_geneid',
            u'K_glpX_uniprotid', u'K_glpX_pos1', u'K_glpX_pos2', u'K_glpX_direct',
            u'predict_O_result', u'predict_K_result',
            u'Predict_O_sero',u'Predict_K_sero',
            u'New_serotype']
    merge_final_output = merge_output.reindex(columns = columns_list)
    merge_final_output.fillna(u'NULL',inplace=True) ##2021-1-18
    #print merge_final_output
    

    #Add_predict_result and if New sero.
      ##A.add_O
        ###ONE judge
    #print("-"*50)
    merge_final_output.loc[merge_final_output[u'O_coaD_contig']==u'NULL',u'Predict_O_sero'] = "One"
    #print merge_final_output

        ###ONT
           ###1.extract but no blast
    merge_final_output.loc[(merge_final_output[u'O_coaD_contig']!=u'NULL')&(merge_final_output[u'predict_O_result']==u'NULL'),u'Predict_O_sero'] = "Ont" #####test is OK!!!
    #print merge_final_output.loc[(merge_final_output[u'O_coaD_contig']!=u'NULL') & (merge_final_output[u'predict_O_result']==u'NULL'),u'Predict_O_sero']
    #print merge_final_output
    #print (1111111111111111111)
          ###2.extract and _a|_b blast.
           ##included in next module


        ##Oxx judge
       #test line
    #merge_final_output[u'predict_O_result'][2]="O10000,O34_a,O45_a"  ###have only a or b,test is OK!!!
    #merge_final_output[u'predict_O_result'][2]="O10000,O34_a,O45_a,O23_b" ##have dif a and b, test is OK!!!
    #merge_final_output[u'predict_O_result'][2]="O10000,O34_a,O45_a,O23_b,O34_b" #have same a and b, test is OK!!
           ###1.contig not NULL & not contain _a|_b----test is OK!!
    merge_final_output.loc[(merge_final_output[u'O_coaD_contig'] != "NULL") & (~merge_final_output[u'predict_O_result'].str.contains("_a|_b")),u'Predict_O_sero'] = merge_final_output.loc[(merge_final_output[u'O_coaD_contig'] != "NULL") & (~merge_final_output[u'predict_O_result'].str.contains("_a|_b")),u'predict_O_result'] 
          
           ###2.congtig not NULL & contain _a|_b & union!=kong--->non _a|_b + union
    row_num =  len(merge_final_output[u'strain_name'])
    
    for i in range(0,row_num):
        i=int(i)
        set_a=[]
        set_b=[]
        set_OK=[]
        if (merge_final_output[u'O_coaD_contig'][i] != "NULL") and ("_a" in merge_final_output[u'predict_O_result'][i] or "_b" in merge_final_output[u'predict_O_result'][i]):#test is OK!!!!
            ele_list = merge_final_output[u'predict_O_result'][i].split(",")
            set_a = [z.split("_")[0] for z in ele_list if "_a" in z]
            set_b = [z.split("_")[0] for z in ele_list if "_b" in z]
            set_OK = [z for z in ele_list if "_a" not in z and "_b" not in z]
           # print set_a
           # print set_b
           # print set_OK
                 #####a and b not two
            if (len(set_a)==0 or len(set_b)==0) and len(set_OK)==0:
                #print "aaa"
                merge_final_output[u'Predict_O_sero'][i]="Ont"
               #merge_final_output.iloc[i,merge_final_output.columns.get_loc(u'Predict_O_sero')]="Ont"
            if (len(set_a)==0 or len(set_b)==0) and len(set_OK)!=0:
                #print type(merge_final_output.columns.get_loc(u'Predict_O_sero'))
                #print type(i)
                #print i
                merge_final_output[u'Predict_O_sero'][i]=",".join(set_OK)
                #merge_final_output.iloc[i,merge_final_output.columns.get_loc(u'Predict_O_sero')]=",".join(set_OK)
                #print "bbb"
                #####a and b not null
            if (len(set_a)!=0 and len(set_b)!=0):
                union_ele = list(set(set_a) & (set(set_b)))
                result_list = union_ele+set_OK
                if len(result_list)==0:
                    merge_final_output[u'Predict_O_sero'][i]="Ont"
                    #print "ccc"
                else:
                    merge_final_output[u'Predict_O_sero'][i]=",".join(result_list)
                    #print "ddd"
                     

    #print(merge_final_output)
            
    # --------------------------------------------------------------------------
      ##B.add_K
        ###KNE judge
    #print("*"*50)
    merge_final_output.loc[merge_final_output[u'K_glpX_contig'] == "NULL",u'Predict_K_sero'] = "Kne"
    #print(merge_final_output)
   
        ###KNT
           ###1.extract but no blast
    merge_final_output.loc[(merge_final_output[u'K_glpX_contig']!=u'NULL')&(merge_final_output[u'predict_K_result']==u'NULL'),u'Predict_K_sero'] = "Knt"
         
           ###2.extract and _a|_b blast.
            ##included in next module
 
        ###Kxx judge
    #merge_final_output[u'predict_K_result'][0]="K10000,K34_a,K45_a"  ###have only a or b,test is OK!!!
    #merge_final_output[u'predict_K_result'][0]="K10000,K34_a,K45_a,K23_b" ##have dif a and b, test is OK!!!
    #merge_final_output[u'predict_K_result'][0]="K10000,K34_a,K45_a,K23_b,K34_b" #have same a and b, test is OK!!
           ###.1.contig not NULL & not contain _a|_b
    merge_final_output.loc[(merge_final_output[u'K_glpX_contig'] != "NULL") & (~merge_final_output[u'predict_K_result'].str.contains("_a|_b")),u'Predict_K_sero'] = merge_final_output.loc[(merge_final_output[u'K_glpX_contig'] != "NULL") & (~merge_final_output[u'predict_K_result'].str.contains("_a|_b")),u'predict_K_result']


           ###.2.contig not NULL & contain _a|_b & union!=kong----->non_a|_b + union
    row_num =  len(merge_final_output[u'strain_name'])

    for i in range(0,row_num):
        i=int(i)
        set_a=[]
        set_b=[]
        set_OK=[]
        if (merge_final_output[u'K_glpX_contig'][i] != "NULL") and ("_a" in merge_final_output[u'predict_K_result'][i] or "_b" in merge_final_output[u'predict_K_result'][i]):
            ele_list = merge_final_output[u'predict_K_result'][i].split(",")
            set_a = [z.split("_")[0] for z in ele_list if "_a" in z]
            set_b = [z.split("_")[0] for z in ele_list if "_b" in z]
            set_OK = [z for z in ele_list if "_a" not in z and "_b" not in z]
            #print set_a
            #print set_b
            #print set_OK
                #####a and b not two
            if (len(set_a)==0 or len(set_b)==0) and len(set_OK)==0:
                #print "aaa"
                merge_final_output[u'Predict_K_sero'][i]="Knt"
            if (len(set_a)==0 or len(set_b)==0) and len(set_OK)!=0:
                #print type(merge_final_output.columns.get_loc(u'Predict_K_sero'))
                #print type(i)
                #print i
                merge_final_output[u'Predict_K_sero'][i]=",".join(set_OK)
                #print "bbb"

                ###a and b not null
            if (len(set_a)!=0 and len(set_b)!=0):
                union_ele = list(set(set_a) & (set(set_b)))
                result_list = union_ele+set_OK
                if len(result_list)==0:
                    merge_final_output[u'Predict_K_sero'][i]="Knt"
                    #print "ccc"
                else:
                    merge_final_output[u'Predict_K_sero'][i]=",".join(result_list)
                    #print "ddd"
    #print(merge_final_output)
      
    ############debug
    merge_final_output.loc[(merge_final_output[u'O_coaD_contig']!=u'NULL')&(merge_final_output[u'predict_O_result']==u'NULL'),u'Predict_O_sero'] = "Ont" #####test is OK!!!
    merge_final_output.loc[(merge_final_output[u'K_glpX_contig']!=u'NULL')&(merge_final_output[u'predict_K_result']==u'NULL'),u'Predict_K_sero'] = "Knt"




    #add New serotype marker
    sero_GB_list = []
    sero_f_dir = program_dir+"source/GB4789-2013_VP_serotype"
    with open (sero_f_dir) as sero_f:
        sero_GB_list = [line.strip() for line in sero_f]
    
    for i in range(0,row_num):
        O_serogroup = merge_final_output[u'Predict_O_sero'][i]
        K_serogroup = merge_final_output[u'Predict_K_sero'][i]
        O_match = re.findall(r'^O[0-9]+',O_serogroup)
        K_match = re.findall(r'^K[0-9]+',K_serogroup)
        #print O_match
        #print K_match
        if len(O_match)>=1 and len(K_match)>=1:
            if O_serogroup+":"+K_serogroup in sero_GB_list:
                merge_final_output[u'New_serotype'][i]="Exist"
            else:
                merge_final_output[u'New_serotype'][i]="New"


    merge_final_output.rename(columns={u'predict_O_result':u'O_Spec_Gene',u'predict_K_result':u'K_Spec_Gene'},inplace=True)
    #colname_list = list(merge_final_output.columns)
    #colname_list
    #print (merge_final_output)
    #print list(merge_final_output.columns)

    #########test##########
    #merge_final_output.iloc[3,27]= "O12"
    #merge_final_output.iloc[1,27]= "O12"
    #merge_final_output.iloc[2,27]= "O7"
    #merge_final_output.iloc[2,28]= "K4"
    #print merge_final_output.iloc[3,27]
    ###################################
    #format the output of predcit_O_sero and predict_K_sero-2011-6-22 
       #method:conditonal select and give value by pandas loc,test is OK!
    with open(program_dir+"source/output_QC.tab") as f_QC:
        O_QC_list = [a.strip() for a in f_QC if a[0]=="O"]
    with open(program_dir+"source/output_QC.tab") as f_QC:
        K_QC_list = [b.strip() for b in f_QC if b[0]=="K"]
    
    for ele in O_QC_list:
        #print(merge_final_output[merge_final_output[u'Predict_O_sero']==ele])
       # merge_final_output[merge_final_output[u'Predict_O_sero']==ele].loc[:,u'Predict_O_sero']="p"+ele
        merge_final_output.loc[merge_final_output[u'Predict_O_sero']==ele,'Predict_O_sero'] = "p"+ele
    
    for ele in K_QC_list:
        merge_final_output.loc[merge_final_output[u'Predict_K_sero']==ele,'Predict_K_sero'] = "p"+ele

    ################################################


    print (merge_final_output) 

    #output to result dir
    dir_cmd = "mkdir -p %s"%(o_dir_all+"04.predict_result")
    os.system (dir_cmd)
    merge_final_output.to_excel(o_dir_all+"04.predict_result/all_strain_predict_result.xlsx")
    
    #print merge_final_output

##############################main program#######################
try:
    if (in_dir != None) and (prokka_dir != None):
        raise RuntimeError
except RuntimeError as e:
    logger.error("Parameter is Error! You can't set the -i and -p parameters at the same time!")
    os._exit(0)

if in_dir != None:
    logger.info("###########################################################################################")
    logger.info("1.genome annote by Prokka begin!")
    genome_annote(thread_num=thread_num)
    logger.info("1.genome annote is OK")
    
    logger.info ("###########################################################################################")
    logger.info("2.blastn to find specific gene begin !")
    query_spec_gene(in_dir=o_dir_all+"01.annote/",in_spec=program_dir+"source/spec_gene_final.fasta")
    logger.info("2.blastn is OK")

    logger.info ("###########################################################################################")
    logger.info ("3.extract gene cluster and border analyse begin !")
    border_analyze(in_dir=o_dir_all+"01.annote/")
    logger.info("3.border analyse is OK")

    logger.info("###########################################################################################")
    logger.info("4.predict O and K serogroup begin !")
    combine_02_03(in_blastn = o_dir_all+"02.blastn/filter_blastn.xlsx",in_O_info = o_dir_all+"03.border_analyze/all_strain_O_info.xlsx", in_K_info = o_dir_all+"03.border_analyze/all_strain_K_info.xlsx")
    logger.info("all is OK")

if prokka_dir != None:
   logger.info("###########################################################################################")
   logger.info("1.copying prokka results to 01.annote/ begin")
   copy_prokka_result()
   logger.info ("1.copying prokka results is OK")
   
   logger.info ("###########################################################################################")
   logger.info("2.blastn to find specific gene begin !")                                                                                                               
   query_spec_gene(in_dir=o_dir_all+"01.annote/",in_spec=program_dir+"source/spec_gene_final.fasta")
   logger.info ("2.blastn is OK")

   logger.info ("###########################################################################################")
   logger.info ("3.extract gene cluster and border analyse begin !")
   border_analyze(in_dir=o_dir_all+"01.annote/")
   logger.info("3.border analyse is OK")
   
   logger.info ("###########################################################################################")
   logger.info ("4.predict O and K serogroup begin !") 
   combine_02_03(in_blastn = o_dir_all+"02.blastn/filter_blastn.xlsx",in_O_info = o_dir_all+"03.border_analyze/all_strain_O_info.xlsx", in_K_info = o_dir_all+"03.border_analyze/all_strain_K_info.xlsx")
   logger.info("all is OK")
