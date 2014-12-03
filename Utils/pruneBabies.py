#! /usr/bin/python

from ROOT import TFile, TTree




def pruneBaby( fname, dir, prunedir, pruneBranches ):

  fullname = dir+"/"+fname
  file = TFile(fullname)
  tree = file.Get("mt2")

  print "-> Pruning " + fullname

  
  for i in pruneBranches:
    tree.SetBranchStatus( i+"_*", 0 )
    tree.SetBranchStatus( "n"+i, 0 )

  newfilename = fname.split(".root")[0] + "_prune.root"
  newfile = TFile(prunedir+"/"+newfilename, "recreate")
  print "newfile: " + newfile.GetName()
  newtree = tree.CloneTree()
  newtree.Write("", 5)
  newfile.Close()



if __name__ == '__main__':

   import os
   import sys

   from optparse import OptionParser

   parser = OptionParser()
   parser.usage = ""
   parser.add_option("-f","--filter", dest="filter",
                      default="",
                      help="prune only selected dataset")

   (options,args) = parser.parse_args()
   if len(args)==0:
     print "ERROR! Script needs at least one argument to run:"
     print "   python pruneBabies.py [input_directory] [output_dir] [configFile=\"pruneBranches.txt\"]"
     exit()

   dir = args[0]
   outdir = args[1]
   pruneBranches = args[2]
   if len(args)>3:
     print "WARNING! Passed more than three arguments to script! Will ignore the additional ones."


   if options.filter !="" :
     print "-> Pruning only files containing: " + str(options.filter)

   pruneBranches = pruneBranches.replace(",", " ")
   pruneBranches = pruneBranches.split()

   print "-> Will prune these branches: " + str(pruneBranches)


   prunedir = outdir
   os.system("mkdir -p " + prunedir)
   print "prunedir: " + prunedir
   print "dir: " + dir

   files = os.listdir(dir)

   for f in files:
     print "found : " + str(f)
     if ".root" in f:
       print "ok got it"
       if options.filter in f:
         pruneBaby(f, dir, prunedir, pruneBranches)

   #print "Find pruned babies in " + prunedir

