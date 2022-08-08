import os
import boto3
import pandas as pd



s3 = boto3.client('s3')

def bucketList():
    response = s3.list_buckets()
    for bucket in response['Buckets']:
        print(f'  {bucket["Name"]}')



#bucketList()

#print(os.getcwd())

# read data from s3
bucket = 'abombin' 
obj = s3.get_object(Bucket= bucket, Key= 'df.csv') 
# get object and file (key) from bucket

df = pd.read_csv(obj['Body']) # 'Body' is a key word
