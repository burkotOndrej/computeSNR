## computeSNR

```
                                 _       _____ _   _ ______ 
                                | |     /  ___| \ | || ___ \
  ___ ___  _ __ ___  _ __  _   _| |_ ___\ `--.|  \| || |_/ /
 / __/ _ \| '_ ` _ \| '_ \| | | | __/ _ \`--. \ . ` ||    / 
| (_| (_) | | | | | | |_) | |_| | ||  __/\__/ / |\  || |\ \ 
 \___\___/|_| |_| |_| .__/ \__,_|\__\___\____/\_| \_/\_| \_|
                    | |                                     
                    |_|    
```
TODO - add some description

Python script used for making 3D ROIs.

These ROIs can be used for computing SNR from MRI image.

### Usage

1. Install

```shell
git clone https://github.com/Razerino/computeSNR.git
```

2. Create `venv` and install requirements

3. Run the script:

```shell
python3 create_noise_ROIs.py -i /home/<user>/sub-001_T1w.nii.gz -visualise 1
```