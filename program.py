#main program: recognize input, use function.
#!/usr/bin/env python
# coding=utf-8

#main program: recognize input, use function.

import argparse
import os
import glob
import time
import re
import pandas as pd
import numpy as np
#parameter input and recognize
parser = argparse.ArgumentParser(description="Prediction for Serotype of Vibrio parahaemolyticus")

parser.add_argument("-i","--in_dir", metavar="<dir>",
                                 required=True,action="store",dest="in_dir",
                                 help="a director that contains all genome assemble fasta file")
parser.add_argument("-o","--out_dir", metavar="<dir>",
                                  required=True,action="store",dest="out_dir",
                                 help="a direcotr that generate analyze result")

parser.add_argument("-n","--thread_num", metavar="<int>",
                                  required=False,action="store",dest="thread_num",                                                              help="set the thread number when genome annotate")

args = parser.parse_args()
in_dir = args.in_dir
o_dir= args.out_dir
thread_num = args.thread_num


#####################################normalize user input
cwd = os.getcwd()

if in_dir[0] !="/":# -i transfer to absolute dir
    in_dir = cwd+"/"+in_dir

if in_dir[-1] != "/":
    in_dir = in_dir+"/"

if o_dir[0] != "/": # -o transfer to absolute dir
    o_dir = cwd+"/" + o_dir

if o_dir[-1] != "/":
    o_dir = o_dir+"/"
os.system("mkdir -p "+o_dir+"serotype_predict")

o_dir_all = o_dir+"serotype_predict/"


if thread_num == None:
    thread_num = 1
else:
    thread_num = int(thread_num)

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
     
    ### make prokka shell in one file
    fna_list = os.listdir(in_dir)
    
    os.system("mkdir -p "+o_dir_all+"middle_file/")

    for i in fna_list:
        prokka_cmd = "echo \"prokka --force --outdir %s --fast --quiet --prefix %s --gcode 11 --locustag %s %s\">>%s"%(annote_dir+i,i,i,in_dir+i,o_dir_all+"middle_file/prokka.sh")
        os.system(prokka_cmd)
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
    print "run cmd:"+split_cmd

    prokka_thread_cmd_file = os.listdir(o_dir_all+"middle_file/")
    for i in prokka_thread_cmd_file:
        os.system("sh %s &"%(i))
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
    print "prokka module is OK!!!"
            
    os.chdir (python_dir)
#############################################################################################
#function2:
    #introduce:input the prokka result dir, specgene fasta file,  and return the blast m8 file.
    #middle progress:1.make blast db
                   # 2.blast
                   # 3.m8 file coverage calculate,filter.first line add.
                   # 4.data toushi sheet.
def query_spec_gene(in_dir=o_dir_all+"01.annote/",in_spec=cwd+"/source/spec_gene_final.fasta"):
    #1.cat all genenome fasta to 1.
    cat_cmd = "cat %s*/*ffn>%s" % (in_dir , o_dir_all+"middle_file/all_genome.fasta")
    os.system (cat_cmd)#OK!! have test!
    os.system ("mkdir -p %s02.blastn"%(o_dir_all))
    print cat_cmd
    
    #2.blast,calculate,filter,data toushi sheet
    
    sh_cmd = "python ./scripts_of/s12_mkdb_blastn.py -idb %s -iq %s -o %s -iL %s" % (o_dir_all+"middle_file/all_genome.fasta", cwd+"/source/spec_gene_final.fasta", o_dir_all+"02.blastn/blastn.m8", cwd+"/source/spec_gene_final.fasta")
    print sh_cmd
    os.system (sh_cmd)
    print "blastn OK!"
    
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
        grep_gff_cmd = "cat %s|grep Prodigal:2.6 |grep %s >%s" % (in_dir+"*/*gff",gene,o_dir_all+"03.border_analyze/"+gene)
        os.system(grep_gff_cmd)
        print grep_gff_cmd
        
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
    hldD = one_border_process(gene="P67911")#the border gene of O and K,hldD
    glpX = one_border_process(gene="P0A9C9")# the border gene of K,glpX
    #test = pd.DataFrame([["05-3133.fna","LRAI01000025.2","05-3133.fna_02345","B7V2S6","174654","175136","+"]])#######test
    #coaD = coaD.append(test,ignore_index=True)########test

    #merge to border table, and rm no on contig. #test! work well!
    coaD_hldD_df = pd.merge(coaD, hldD, left_on=[0],right_on=[0],how="outer")
    coaD_hldD_df = coaD_hldD_df[coaD_hldD_df["1_x"] == coaD_hldD_df["1_y"] ]
    print coaD_hldD_df

    hldD_glpX_df = pd.merge(hldD, glpX, left_on=[0],right_on=[0],how="outer")
    hldD_glpX_df = hldD_glpX_df[hldD_glpX_df["1_x"] == hldD_glpX_df["1_y"]]
    print hldD_glpX_df

    #left merge the whole strain number, and generate the final table.
    all_strain_list = os.listdir(in_dir)
    all_strain_df = pd.DataFrame(all_strain_list)

    all_strain_O_df = pd.merge(all_strain_df, coaD_hldD_df, left_on=[0],right_on=[0],how="left")
    print all_strain_O_df
    all_strain_K_df = pd.merge(all_strain_df, hldD_glpX_df, left_on=[0],right_on=[0],how="left")
    print all_strain_K_df

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
    print filter_blastn_df

    path_2 = unicode(in_O_info,"utf8")
    O_cluter_df = pd.read_excel(path_2, sheet_name = "Sheet1")#have test, OK!
    O_cluter_df.columns = ["O_strain_id","O_coaD_contig","O_coaD_geneid","O_coaD_uniprotid","O_coaD_pos1","O_coaD_pos2","O_coaD_direct","O_hldD_contig","O_hldD_geneid","O_hldD_uniprotid","O_hldD_pos1","O_hldD_pos2","O_hldD_direct"]
    print O_cluter_df

    path_3 = unicode(in_K_info,"utf8")
    K_cluter_df = pd.read_excel(path_3, sheet_name = "Sheet1")#have test, OK!
    K_cluter_df.columns = ["K_strain_id","K_hldD_contig","K_hldD_geneid","K_hldD_uniprotid","K_hldD_pos1","K_hldD_pos2","K_hldD_direct","K_glpX_contig","K_glpX_geneid","K_glpX_uniprotid","K_glpX_pos1","K_glpX_pos2","K_glpX_direct"]
    print K_cluter_df

    #2.filter blastn result(in K cluster or O cluster)
    filter_blastn_O_cluster = pd.merge(filter_blastn_df, O_cluter_df, left_on=["subject_strain_id"],right_on=["O_strain_id"],how="left")
    filter_blastn_OK_cluster = pd.merge(filter_blastn_O_cluster, K_cluter_df, left_on=["subject_strain_id"],right_on=["K_strain_id"],how="left")
    print filter_blastn_OK_cluster
    #filter_blastn_OK_cluster.to_excel(o_dir_all+"03.border_analyze/filter_blastn_OK_cluster.xlsx")
    
    #2.A O predict result:
    blastn_O_df = filter_blastn_OK_cluster[filter_blastn_OK_cluster["query_id"].apply(lambda x:x[0])=="O"]
    blastn_O_dropna_df = blastn_O_df.dropna(axis=0,subset=["O_coaD_geneid"])#rm no O cluster OK!!!
    print blastn_O_dropna_df
    blastn_O_dropna_df["subject_locus"] = blastn_O_dropna_df["subject_id"].apply(lambda x:int(x.split("_")[-1]))
    blastn_O_dropna_df["O_coaD_locus"] = blastn_O_dropna_df["O_coaD_geneid"].apply(lambda x:int(x.split("_")[-1]))
    blastn_O_dropna_df["O_hldD_locus"] = blastn_O_dropna_df["O_hldD_geneid"].apply(lambda x:int(x.split("_")[-1]))
    blastn_O_dropna_df["O_max"] = blastn_O_dropna_df[["O_coaD_locus","O_hldD_locus"]].max(axis = 1)
    blastn_O_dropna_df["O_min"] = blastn_O_dropna_df[["O_coaD_locus","O_hldD_locus"]].min(axis = 1)
    #print blastn_O_dropna_df
    blastn_O_dropna_df = blastn_O_dropna_df[(blastn_O_dropna_df["subject_locus"]<=blastn_O_dropna_df["O_max"]) & (blastn_O_dropna_df["subject_locus"]>=blastn_O_dropna_df["O_min"])]#rm the subject not in O cluster, have test, work well!!!
    print blastn_O_dropna_df
    #blastn_O_dropna_df.to_excel(o_dir_all+"03.border_analyze/blastn_O_dropna_df.xlsx")

    #2.B K predict result:
    blastn_K_df = filter_blastn_OK_cluster[filter_blastn_OK_cluster["query_id"].apply(lambda x:x[0])=="K"]
    blastn_K_dropna_df = blastn_K_df.dropna(axis=0,subset=["K_glpX_geneid"])#rm no K cluster OK!!!
    print blastn_K_dropna_df
    blastn_K_dropna_df["subject_locus"] = blastn_K_dropna_df["subject_id"].apply(lambda x:int(x.split("_")[-1]))
    blastn_K_dropna_df["K_hldD_locus"] = blastn_K_dropna_df["K_hldD_geneid"].apply(lambda x:int(x.split("_")[-1]))
    blastn_K_dropna_df["K_glpX_locus"] = blastn_K_dropna_df["K_glpX_geneid"].apply(lambda x:int(x.split("_")[-1]))
    blastn_K_dropna_df["K_max"] = blastn_K_dropna_df[["K_hldD_locus","K_glpX_locus"]].max(axis = 1)
    blastn_K_dropna_df["K_min"] = blastn_K_dropna_df[["K_hldD_locus","K_glpX_locus"]].min(axis = 1)
    #print blastn_K_dropna_df

    blastn_K_dropna_df = blastn_K_dropna_df[(blastn_K_dropna_df["subject_locus"]<=blastn_K_dropna_df["K_max"]) & (blastn_K_dropna_df["subject_locus"]>=blastn_K_dropna_df["K_min"])]#rm the subject not in K cluster ,have test, work well!!!!
    print blastn_K_dropna_df

    #blastn_K_dropna_df.to_excel(o_dir_all+"03.border_analyze/blastn_K_dropna_df.xlsx")

     #3 long to pivot
       #3.A O
    blastn_O_pivot_df = blastn_O_dropna_df.pivot_table (index = "subject_strain_id", columns = "query_id", values = "subject_id", aggfunc = "count")    
    print blastn_O_pivot_df
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
    print final_predict_list
    O_predict_df = pd.DataFrame(final_predict_list,columns = ["strain_id","predict_O_result"])
    print O_predict_df


       #3.B K
    blastn_K_pivot_df = blastn_K_dropna_df.pivot_table (index = "subject_strain_id", columns = "query_id", values = "subject_id", aggfunc = "count")
    print blastn_K_pivot_df
    row_list2 = [d for d in blastn_K_pivot_df.index]
    row_num2 = len(row_list2)
    final_predict_list2 = []
    for e in range(0,row_num2):
        line2 = blastn_K_pivot_df.iloc[e]
        line2_dropna = line2.dropna()
        predict_list = [f.split("*")[0] for f in line2_dropna.index]
        predict_result = ",".join(predict_list)
        print predict_result
        out_strain_name = row_list2[e]
        final_predict_list2.append([out_strain_name,predict_result])
    print final_predict_list2
    K_predict_df = pd.DataFrame(final_predict_list2,columns = ["strain_id","predict_K_result"])
    print K_predict_df

    #4.final output result
    all_strain_list = os.listdir(o_dir_all+"/01.annote")    
    print all_strain_list
    all_stain_df = pd.DataFrame([[ele] for ele in all_strain_list],columns=["strain_name"])
    print K_predict_df
    print "-"*50
    merge_df1 = pd.merge(all_stain_df, O_cluter_df, left_on=["strain_name"],right_on=["O_strain_id"],how="left")
    merge_df2 = pd.merge(merge_df1,K_cluter_df,left_on=["strain_name"],right_on=["K_strain_id"],how="left")
    print merge_df2
    
    merge_df3 =  pd.merge(merge_df2,O_predict_df,left_on=["strain_name"],right_on=["strain_id"],how="left")
    print merge_df3

    merge_output = pd.merge(merge_df3,K_predict_df,left_on=["strain_name"],right_on=["strain_id"],how="left")
    print merge_output
    print merge_output.columns
    columns_list = [u'strain_name', u'O_coaD_contig', u'O_coaD_geneid',u'O_coaD_uniprotid', u'O_coaD_pos1', u'O_coaD_pos2', u'O_coaD_direct',
            u'O_hldD_contig', u'O_hldD_geneid', u'O_hldD_uniprotid', u'O_hldD_pos1',
            u'O_hldD_pos2', u'O_hldD_direct', u'K_hldD_contig',
            u'K_hldD_geneid', u'K_hldD_uniprotid', u'K_hldD_pos1', u'K_hldD_pos2',
            u'K_hldD_direct', u'K_glpX_contig', u'K_glpX_geneid',
            u'K_glpX_uniprotid', u'K_glpX_pos1', u'K_glpX_pos2', u'K_glpX_direct',
            u'predict_O_result', 
            u'predict_K_result']
    merge_final_output = merge_output.reindex(columns = columns_list)
    print merge_final_output

    #output to result dir
    dir_cmd = "mkdir -p %s"%(o_dir_all+"04.predict_result")
    os.system (dir_cmd)
    merge_final_output.to_excel(o_dir_all+"04.predict_result/all_strain_predict_result.xlsx")
##############################main program#######################
#genome_annote(thread_num=thread_num)
#query_spec_gene(in_dir=o_dir_all+"01.annote/",in_spec=cwd+"/source/spec_gene_final.fasta")
#border_analyze(in_dir=o_dir_all+"01.annote/")
combine_02_03(in_blastn = o_dir_all+"02.blastn/filter_blastn.xlsx",in_O_info = o_dir_all+"03.border_analyze/all_strain_O_info.xlsx", in_K_info = o_dir_all+"03.border_analyze/all_strain_K_info.xlsx")
