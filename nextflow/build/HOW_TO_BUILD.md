# How to Build Apptainer Containers on Windows Desktop

## Section 1: One-Time Setup (do this once ever)

### Install WSL with Ubuntu
Press the Windows key, search "PowerShell", right-click it and select "Run as Administrator". Then run:
```
wsl --install -d Ubuntu
```
Reboot when prompted.

### Install Apptainer in WSL
Press the Windows key, search "Ubuntu", click to open it. Then run:
```
sudo apt update && sudo apt upgrade -y
sudo apt install -y software-properties-common
sudo add-apt-repository -y ppa:apptainer/ppa
sudo apt update
sudo apt install -y apptainer
```

---

## Section 2: Before Each Build

### Navigate to your build folder
Press the Windows key, search "Ubuntu", click to open it. Then run:
```
cd "/mnt/c/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/build"
```

### What files must be in the build folder

| Container        | Required file                                      | Where to get it                                                                 |
|------------------|----------------------------------------------------|---------------------------------------------------------------------------------|
| cellranger       | cellranger-10.0.0.tar.gz                           | https://www.10xgenomics.com/support/software/cell-ranger/downloads              |
| spaceranger      | spaceranger-4.1.0.tar.gz                           | https://www.10xgenomics.com/support/software/space-ranger/downloads             |
| gdc_client       | gdc-client_2.3_Ubuntu_x64-py3.8-ubuntu-20.04.zip   | https://gdc.cancer.gov/access-data/gdc-data-transfer-tool                      |
| omics            | omics_r4.4.3.def, omics_cache_omnipath_data.R      | —                                                                               |

> **Note:** 10x Genomics download URLs are signed and expire after a few days.
> If download fails with a 403 error, go to the link above, right-click the download
> button, copy the link address, and use that fresh URL.

### Verify files are present before building
```
ls -lh
```

---

## Section 3: Build Commands

Always use `2>&1 | tee build.log` so output is shown on screen AND saved to a log file.
If the build fails, open `build.log` and Ctrl+F for `ERROR` to find the exact cause.

### Cell Ranger
```
sudo apptainer build cellranger_10.0.0.sif cellranger.def 2>&1 | tee cellranger.log
```

### Space Ranger
```
sudo apptainer build spaceranger_4.1.0.sif spaceranger.def 2>&1 | tee spaceranger.log
```

### GDC Client
```
sudo apptainer build gdc_client.sif gdc_client.def 2>&1 | tee gdc_client.log
```

### Omics
First, run omics_cache_omnipath_data.R. Then, run below command
```
sudo apptainer build omics_r4.4.3.sif omics_r4.4.3.def 2>&1 | tee omics.log
sudo apptainer build omics_py3.12.sif omics_py3.12.def 2>&1 | tee build.log
```

Typical build times:
- cellranger / spaceranger: 15-30 minutes
- gdc_client: 5 minutes
- omics: 60-120 minutes

---

## Section 4: Upload to HPC

After a successful build, copy the .sif file to HPC:
```
scp cellranger_10.0.0.sif username@hpc.edu:/path/to/containers/
scp spaceranger_4.1.0.sif username@hpc.edu:/path/to/containers/
scp gdc_client.sif username@hpc.edu:/path/to/containers/
scp omics_2026.04.29.sif username@hpc.edu:/path/to/containers/
```

---

## Section 5: Troubleshooting

**"Permission denied"**
Always use `sudo` before apptainer commands.

**"No space left on device"**
WSL has a default disk limit. To increase it, create or edit `C:\Users\kailasamms\.wslconfig`:
```
[wsl2]
diskSize=100GB
```
Then restart WSL: `wsl --shutdown` in PowerShell, then reopen WSL.

**"Can't find files"**
Access Windows files from WSL using `/mnt/c/` prefix, not `C:\`.

**Build failed — how to find the error**
```
grep -i "error" build.log
```
Or open `build.log` in any text editor and Ctrl+F for `ERROR`.

**10x download URL expired (403 Forbidden)**
Go to the 10x downloads page, right-click the download button, copy the fresh signed URL, and download manually:
```
wget -c -O cellranger-10.0.0.tar.gz "PASTE_FRESH_URL_HERE"
```

**"command not found" after build**
Run the %test block manually to verify:
```
sudo apptainer test cellranger_10.0.0.sif
```
