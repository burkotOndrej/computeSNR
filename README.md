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
Python script used for making 2D or 3D background noise ROIs and SNR calculation.

Optional arguments:

    -shiftx : Shift of the ROIs in X axis from the image borders.
    -shifty : Shift of the ROIs in Y axis from the image borders.
    -shiftunits : Shift units in pixels (pix) or percentage (per).
    -size : Size of ROIs in X and Y axis in pixels.
    -visualise : Visualisation of created ROI. 0 - do not visualise, 1 - visualise'
    -outpath : Path for saving nii file with created ROIs.

### Usage

1. Install

```shell
$ git clone https://github.com/Razerino/computeSNR.git
```

2. Navigate to cloned folder with script

```shell
$ cd ~/<usr>/<directory>/
```

3. Create `venv` and install requirements

```shell
$ python3 -m venv venv

$ source venv/bin/activate

$ pip install -r requirements.txt
```

3. Run the script:

```shell
$ python3 computeSNR.py -i /home/<user>/sub-001_T1w.nii.gz -m /home/<user>/sub-001_T1w_roi_mask.nii.gz 
```
Note: It is always good practice to manually input also optional arguments (at least `-shiftx` and `-shifty`)
because then you will be sure, that noise mask is generated on right places.

4. You can show help for script

```shell
$ computeSNR.py -h
```