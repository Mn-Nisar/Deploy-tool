import json
import requests
import boto3
import pandas as pd
import io
import os 

S3_BUCKET=''
AWS_ACCESS_KEY_ID=''
S3_SECRET_KEY=''

mybucket = S3_BUCKET

s3 = boto3.resource('s3',
         aws_access_key_id=AWS_ACCESS_KEY_ID,
         aws_secret_access_key= S3_SECRET_KEY)

def create_s3_object(df,object_id):

    csv_buffer = io.StringIO()
    df.to_csv(csv_buffer, index = False)
    
    try:
        s3.Object(mybucket, object_id).put(Body=csv_buffer.getvalue())
    except:
        return Exception

def batch_correct_limma(df, batches_list, job_id):
    object_id = str(job_id)+"for_batch.csv"
    batches = json.dumps(batches_list)
    create_s3_object(df,object_id)
    params = {'batches':batches, 's3Object':object_id}
    batch_correct_api = requests.post("http://limma_api:8000/remove_batch_effect", params=params)    
    json_data = json.load(io.BytesIO(batch_correct_api.content))
    btc_df = pd.DataFrame(json_data,columns=list(df.columns))
    # delete_s3_object(object_id)
    return btc_df
  

def limma_diff_API(df, matrix,job_id, rename_list):
    object_id = str(job_id)+"for_diffretial.csv"
    print(object_id)
    matrix_li = json.dumps(matrix)
    create_s3_object(df,object_id)
    params = {'matrix_li':matrix_li, 's3Object_for_diff':object_id}

    limma_diff = requests.post("http://limma_api:8000/limma_diff_calc", params=params)    

    json_data = json.load(io.BytesIO(limma_diff.content))
    
    limma_df = pd.DataFrame(json_data)

    limma_df.columns = rename_list

    return limma_df


def getplot_example_s3(plot_type):
    session = boto3.Session(
    aws_access_key_id = AWS_ACCESS_KEY_ID,
    aws_secret_access_key=S3_SECRET_KEY)
    s3 = session.client('s3')

    obj = s3.get_object(Bucket='plotexampleom', Key=plot_type+'.csv')

    df = pd.read_csv(obj['Body'],encoding_errors='ignore')
    return df

