# does not work and logic is unclear without clear product definition/manual revisioning 
def addCodonNumb(annotDf):
  editDf = pd.DataFrame()
  productsList = list(pd.unique(annotDf["All_Products"]))
  for i in productsList:
    print(i)
    print((annotDf.index.to_series().diff().fillna(1) == 1).all())
    curProd = annotDf[annotDf['All_Products'] == i].reset_index().drop(['index'], axis =1)
    prodUn = list(pd.unique(curProd["All_Products"]))
    if (prodUn[0] == "NaN") or (prodUn[0] == "NCR"):
      codList = ["NCR"] * len(curProd)
    else:
      posList = list(range(1, 1 + int((len(curProd.index)/3))))
      codList = []
      for item in posList:
        codList.extend([item]*3)
        
    curProd["Codon#"] = codList
    editDf = pd.concat([editDf, curProd], ignore_index = True)
    
  sortDf = combDf.sort_values(by=['NT']).reset_index().drop(['index'], axis = 1)
  return(sortDf)

annotDf= addCodonNumb(annotDf)
# another approach to try is subset datatable for not NCR or NaN in All_products, count codons, join and reorder


