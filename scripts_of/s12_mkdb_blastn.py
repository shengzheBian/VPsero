#!/usr/bin/env python
# coding=utf-8
#功能：接受一个序列文件，自动对齐进行blastn建库,
      #同时接受一个搜索基因，在建立的数据库中进行搜索。
      #输出一个m8文件，对其进行

import argparse
import os
import pandas as pd
import numpy as np
from Bio import SeqIO
#构建参数输入
parser = argparse.ArgumentParser(description="功能：输入一个fasta文件，自动对其进行blastn建库,根据提供的ref序列blast，并统计其长度。最好在建库目录下运行。问题为python统计结果输出不出来。命令示例：python xxx.py -idb xxx -iq yyy -o 03.blastn/m8/E249_O111.m8")

parser.add_argument("-idb","--in_db", metavar="<file>",
                   required=True,action="store",dest="db_file",
                   help="输入需要建库的fasta文件绝对路径")
parser.add_argument("-iq","--in_query", metavar="<file>",
                    required=True,action="store",dest="query_file",
                   help="输入需要搜索的参考基因序列文件绝对路径")
parser.add_argument("-o","--out_filename", metavar="<file>",
                    required=True,action="store",dest="out_file",
                    help="指定最终生成的m8文件的文件名和地址")

parser.add_argument("-iL","--in_length", metavar="<file>",
                    required=True,action="store",dest="L_file",
                    help="指定需要统计序列长度的fasta文件")
args = parser.parse_args()
db_file= args.db_file
q_file= args.query_file
o_file= args.out_file
L_file=args.L_file

#脚本运行命令日志
mkdb_cmd="formatdb -i "+db_file+" -p F"
q_cmd= "/share/app/blast-2.2.26/bin/blastall -p blastn -e 1e-5 -i "+q_file+" -o "+o_file+" -F F -d "+db_file+" -m 8 -W 11"
qLen_cmd= "python /hwfssz5/ST_INFECTION/Salmonella/liqiwen/bianshengzhe/liliqiang/bianshengzhe/bin/s12_count_fasta_length_biopython.py -i "+L_file
print "将运行的命令如下"
print mkdb_cmd
print q_cmd
print qLen_cmd

#对接收的文件进行建库
os.system(mkdb_cmd)

#搜索需要查找的序列。
os.system(q_cmd)

#去掉工作目录生成的建库日志。
os.system("rm ./formatdb.log")
#统计参考序列长度
print "序列长度统计结果如下："
os.system(qLen_cmd)

#处理m8数据，加一列列名，加一列coverage，filter，数据透视

def m8_to_matrix(in_m8 = o_file, in_length_fasta = L_file):
  #在m8文件中添加列名。
  colnames = "query_id\tsubject_id\tidentity\talignment_length\tmismatched\tgap_opening\tq.start\tq.end\ts.start\ts.end\te-value\tbit_score"
  add_cmd = "sed -i -e '1i %s' %s" % (colnames,in_m8)
  os.system (add_cmd)
  
  #################在m8文件中添加一列qlength和coverage列，过滤。
  #1.读入m8文件。
  m8_df = pd.read_csv(in_m8, sep = '\t')#读入m8文件表,检验了该表格，与m8文件一致无误！
  #m8_df.columns = colnames.split("\t")
  print m8_df.columns
  print m8_df

  #2.生成qlength表，添加到m8表中,并且添加一列coverage。
  length_list = []
  for seq_record in SeqIO.parse(in_length_fasta,"fasta"):
    length_list.append([seq_record.id,len(seq_record)])
    #print seq_record.id,len(seq_record)
  length_df = pd.DataFrame(length_list)#获得refgene 长度df表。
  length_df.columns = ["genename","ref_length"]
  print length_df

  merge_df = pd.merge(m8_df, length_df, left_on=["query_id"],right_on=["genename"],how="left")
  merge_df["coverage"] = merge_df["alignment_length"]/merge_df["ref_length"]
  print merge_df
  
 #3.过滤原始的m8表。
  filter_df = merge_df[(merge_df["identity"]>=90)&(merge_df["coverage"]>=0.8)]#测试工作正常
  print filter_df

 #4.生成数据透视表。
  pivot_df = filter_df.pivot_table(index="subject_id",columns="query_id",values="bit_score",aggfunc="count") #已本地excel测试，该透视表无误一致。 
  print pivot_df

 #输出blastn该步的所有结果文件
  o_dir = "/".join(in_m8.split("/")[0:-1])
  merge_df.to_excel(o_dir+"/raw_blastn.xlsx")#输出未filter前的加了coverage的m8文件。
  filter_df.to_excel(o_dir+"/filter_blastn.xlsx")  
  pivot_df.to_excel(o_dir+"/pivot.xlsx")
 
 #清空middle_file中的中间结果。
  os.system("rm %s" %(o_dir+"/../middle_file/all_genome.fasta*"))
#程序结束

#调用函数
m8_to_matrix(in_m8 = o_file, in_length_fasta = L_file)
print "The blastn analyse is OK!!!"
