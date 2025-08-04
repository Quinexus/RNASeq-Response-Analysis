#!/usr/bin/env python3

import pandas as pd
import os
import subprocess
import ftplib
import sys

ftp_server = 'ftp.ncbi.nlm.nih.gov'
accession_template = "/geo/series/{removed}nnn/{full}/"

def gunzip_file(file):
    subprocess.run(["gunzip", file])

def download_ftp(ftp_dir, ftp, gse, required_ftp_dir):
    ftp.cwd(required_ftp_dir + ftp_dir + "/")
    output_dir = gse + "/" + ftp_dir
    os.makedirs(output_dir, exist_ok=True)

    filenames = ftp.nlst()

    for filename in filenames:
        local_path = os.path.join(output_dir, filename)
        with open(local_path, 'wb') as f:
            ftp.retrbinary(f'RETR {filename}', f.write)
        
        if ".gz" in filename:
            gunzip_file(local_path)

def process_gse(gse):
    os.makedirs(gse, exist_ok=True)
    required_ftp_dir = accession_template.format(removed=gse[:-3], full=gse)

    with ftplib.FTP(ftp_server) as ftp:
        ftp.login()
        ftp.cwd(required_ftp_dir)
        dirs_list = ftp.nlst()
        
        for directory in dirs_list:
            download_ftp(directory, ftp, gse, required_ftp_dir)
    
    return("Processed " + gse)


def main():
    # Check if enough arguments were passed
    if len(sys.argv) < 2:
        print("enter GSE")
        sys.exit(1)
    
    required_gse = sys.argv[1]
    process_gse(required_gse)

    


if __name__ == "__main__":
    main()












