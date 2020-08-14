#!/usr/bin/env python
# coding=utf-8
#功能：输入一个fasta文件head行中的非冗余关键词，返回该基因整个序列。
import argparse

parser = argparse.ArgumentParser(description = "功能：输入一个fasta文件head行中的非冗余关键词文件，返回该基因整个序列。")
parser.add_argument("-q", "--query", metavar="<file>",
                       required=True, action="store",dest="locus_gene_file",
                       help="输入一个包含各搜寻关键词的列表，每个词一行,必须是唯一标识")
parser.add_argument("-db","--database",metavar="<.fasta>",
                     required=True, action="store",dest="fasta_file",
                     help="输入一个用于查找提取序列的fasta格式数据库文件,副溶40K的1420orf序列文件为：/hwfssz5/ST_INFECTION/Salmonella/liqiwen/bianshengzhe/liliqiang/bianshengzhe/VP1pipeexam/key_table_pic_and_control_version/third_serotype/2.41K_info/2.40K_fasta/40K_1420orf_ffn.fasta")

args = parser.parse_args()
query = args.locus_gene_file
db = args.fasta_file

#print query
#print db

##1.1存储传入的基因locus id为列表
locus_gene_file = open(query)
locus_id_list = [i.strip() for i in locus_gene_file]

##1.2读入fasta序列数据库
fasta_file = open(db)
fasta_database = fasta_file.readlines()

#for i in fasta_database:
#    b = i.strip()
#    print b      #测试成功，读取数据库，去掉换行符打印。

#查找并打印一个locus_id
def find_a_locus_id_fasta(locus_id):#####千年大坑！！！for迭代读文件，只能在一个程序中进行一次，写在循环中只能一次有效。
    status = 0
    for line in fasta_database:
        
        if status == 1:
            if line[0] == ">":
                status= 0
            else:
                line_body=line.strip()
                print line_body
        if locus_id in line:
            line_head = line.strip()
            print line_head
            status = 1



for locus_id in locus_id_list:
    find_a_locus_id_fasta(locus_id)
