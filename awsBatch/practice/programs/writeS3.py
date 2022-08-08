import os
import boto3
import pandas as pd
from io import StringIO



s3 = boto3.client('s3')
s3_resource = boto3.resource('s3')


def bucketList():
    response = s3.list_buckets()
    for bucket in response['Buckets']:
        print(f'  {bucket["Name"]}')



#bucketList()

#print(os.getcwd())

# read data from the server
myDf=pd.read_csv('/home/ubuntu/github/shinyApp/shinyVizApps/uploadDatCov/data/metaData.tsv', sep='\t')

#print(df.columns)

# write data to s3
myBucket = 'abombin' 

def writeToS3(bucket, df):
    csv_buffer = StringIO()
    df.to_csv(csv_buffer, index=False)
    s3_resource.Object(bucket, 'metaDat1.csv').put(Body=csv_buffer.getvalue())

bucketList()

writeToS3(myBucket, myDf)


