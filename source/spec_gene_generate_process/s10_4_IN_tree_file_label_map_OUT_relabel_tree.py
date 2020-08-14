#!/usr/bin/env python
# coding=utf-8

import argparse
import os
import re
#功能：输入一个只包含tree文件的目录，不加“/”,以及输入一个label map关系表，原来的标签必须非冗余，即不能一对多映射。然后在当前目录输出所有relabel后的tree文件。
#获取树目录下所有文件名
parser=argparse.ArgumentParser(description = "功能：批量修改nwk文件的label注释。输入一个包含tree文件的目录，和一个tree中老label和要替换的label的对应两列表。在当前目录输出relabel后的各个tree文件")
parser.add_argument("-d","--dir",metavar="<str>",
                    required=True,action="store",dest="input_dir",
                    help="一个只包含要修改的tree文件的目录，必须是绝对路径，最后面没有“/”")
parser.add_argument("-t","--table",metavar="<file>",
                    required=True,action="store",dest="map_dir",
                    help="一个两列表，第一列是树文件中原始的label名称，第二列是想修改为的label名称，中间用tab隔开。前后一一对应，且第一列不能一对多，可以直接由excel复制过来")

args = parser.parse_args()
input_dir = args.input_dir
map_dir = args.map_dir

#input_dir = "/hwfssz5/ST_INFECTION/Salmonella/liqiwen/bianshengzhe/liliqiang/bianshengzhe/VP1pipeexam/10.41K_analyse/s4.orthofinder_result_analyse/s1.over_group_gene_tree_change_label/Recon_Gene_Trees"

tree_file_list = os.listdir(input_dir)

#读入raw_label和new_label映射表为map文件
#map_dir = "/hwfssz5/ST_INFECTION/Salmonella/liqiwen/bianshengzhe/liliqiang/bianshengzhe/VP1pipeexam/10.41K_analyse/s4.orthofinder_result_analyse/s1.over_group_gene_tree_change_label/raw_treelabel_new.txt"
with open(map_dir,"r") as map_dir_f:
    map_dict = {line.strip().split("%")[0]:line.strip().split("%")[1] for line in map_dir_f}
    print map_dict
#ralabel算法核心：1.粗暴法，直接str.replace方法，先搜寻是否有一样的label名字，如果有匹配，就替换为对应的新name。
    
def relabel_tree(file_dir):
    with open(file_dir,"r") as raw_tree_f:#可以正确读取原来树文件。
        raw_tree_str = raw_tree_f.read()
        new_tree_str=raw_tree_str
    for label in map_dict.keys():
        new_tree_str=new_tree_str.replace(label,map_dict[label])
    print new_tree_str
    out_file_name= file_dir.split("/") [-1].split(".")[0]+"_relabel.txt"
    with open(out_file_name,"w") as f:
        print >> f,new_tree_str


for i in tree_file_list:
    file_dir = input_dir+"/"+i
    relabel_tree(file_dir)

