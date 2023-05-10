import pandas as pd
import os
import numpy as np
import pytest

os.chdir('C:/Users/abomb/OneDrive - Emory University/Variant-calling-pipeline/Original_output_files_02032023')

refRes = pd.read_csv("test_results.csv")

resDf = pd.read_csv("filtered.csv")


def test_Type(ref, res, i, testPos):
  # type test
  refType = ref.loc[i, "Type"].lower()
  resType = res.loc[res["POSITION"] == testPos, "Type"].iloc[0].lower()
  try:
    assert refType == resType
  except:
    print("testType" + str(i+1) + " Row  " + str(testPos) + " Fail")
    print(refVar + " " + resVar)

def test_CorrecPos(ref, res, i, testPos):
  refType = ref.loc[i, "Position_corrected"]
  resType = res.loc[res["POSITION"] == testPos, "Position_corrected"].iloc[0]
  try:
    assert refType == resType
  except:
    print("testCorrecPos" + str(i+1) + " Row  " + str(testPos) + " Fail")
    print(refVar + " " + resVar)

def test_VarAlCor(ref, res, i, testPos):
  refVar = ref.loc[i, "VAR-allele_corrected"].lower()
  resVar = res.loc[res["POSITION"] == testPos, "VAR-allele_corrected"].iloc[0].lower()
  try:
    assert refVar== resVar
  except:
    print("VAR-allele_corrected " + str(i+1) + " Row  " + str(testPos) + " Fail")
    print(refVar + " " + resVar)

def test_RefAlCor(ref, res, i, testPos):
  refVar = ref.loc[i, "REF-allele_corrected"].lower()
  resVar = res.loc[res["POSITION"] == testPos, "REF-allele_corrected"].iloc[0].lower()
  try:
    assert refVar== resVar
  except:
    print("REF-allele_corrected " + str(i+1) + " Row  " + str(testPos) + " Fail")
    print(refVar + " " + resVar)

def test_Level(ref, res, i, testPos):
  refVar = ref.loc[i, "Level"].lower()
  resVar = res.loc[res["POSITION"] == testPos, "Level"].iloc[0].lower()
  try:
    assert refVar== resVar
  except:
    print("Level " + str(i+1) + " Row  " + str(testPos) + " Fail")
    print(refVar + " " + resVar)

def test_PosTest(ref, res, i, testPos):
  refVar = ref.loc[i, "Position_test"].lower()
  resVar = res.loc[res["POSITION"] == testPos, "Position_test"].iloc[0].lower()
  try:
    assert refVar== resVar
  except:
    print("Position_test " + str(i+1) + " Row  " + str(testPos) + " Fail")
    print(refVar + " " + resVar)
    
def test_NucleotideChange(ref, res, i, testPos):
  refVar = ref.loc[i, "Nucleotide_Change"].lower()
  resVar = res.loc[res["POSITION"] == testPos, "Nucleotide_Change"].iloc[0].lower()
  try:
    assert refVar== resVar
  except:
    print("Nucleotide_Change " + str(i+1) + " Row  " + str(testPos) + " Fail")
    print(refVar + " " + resVar)

def test_Freq_adj(ref, res, i, testPos):
  refVar = ref.loc[i, "Freq_adj"]
  resVar = res.loc[res["POSITION"] == testPos, "Freq_adj"].iloc[0]
  try:
    assert refVar== resVar
  except:
    print("Freq_adj " + str(i+1) + " Row  " + str(testPos) + " Fail")
    print(refVar + " " + resVar)
    
def test_Region(ref, res, i, testPos):
  refVar = ref.loc[i, "Region"].lower()
  resVar = res.loc[res["POSITION"] == testPos, "Region"].iloc[0].lower()
  try:
    assert refVar== resVar
  except:
    print("Region " + str(i+1) + " Row  " + str(testPos) + " Fail")
    print(refVar + " " + resVar)

def test_Mutation_type(ref, res, i, testPos):
  refVar = ref.loc[i, "TYPE"].lower()
  resVar = res.loc[res["POSITION"] == testPos, "Mutation_type"].iloc[0].lower()
  try:
    assert refVar== resVar
  except:
    print("Mutation_type " + str(i+1) + " Row  " + str(testPos) + " Fail")
    print(refVar + " " + resVar)

def test_AA_change(ref, res, i, testPos):
  refVar = ref.loc[i, "AA_CHANGE"].lower()
  resVar = res.loc[res["POSITION"] == testPos, "AA_change"].iloc[0].lower()
  try:
    assert refVar == resVar
  except:
    print("AA_change " + str(i+1) + " Row  " + str(testPos) + " Fail")
    print(refVar + " " + resVar)

# run all checks
def test_checkMatch(ref, res):
  for i in range(len(ref.index)):
    testPos = ref.loc[i, "POSITION"]
    test_Type(ref, res, i, testPos)
    test_CorrecPos(ref, res, i, testPos)
    test_RefAlCor(ref, res, i, testPos)
    test_VarAlCor(ref, res, i, testPos)
    test_Level(ref, res, i, testPos)
    test_PosTest(ref, res, i, testPos)
    test_NucleotideChange(ref, res, i, testPos)
    test_Freq_adj(ref, res, i, testPos)
    test_Region(ref, res, i, testPos)
    test_Mutation_type(ref, res, i, testPos)
    test_AA_change(ref, res, i, testPos)

test_checkMatch(ref= refRes, res = resDf)

    

