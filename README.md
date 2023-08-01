This is a tiny plot tool of matlab for ulog form PX4, which inspired by [matulog](https://github.com/CarlOlsson/matulog). You need to plot what you want by your shelf.

it have been test on MacOS Monterey 12.3 and win10, matlab 2021b.

Requirements:

1.matlab2016a (or a newer version)

2.python.

3.pyulog is needed for `ulog_matlab_plot` to run, install the latter using ```pip install pyulog``` or ```pip3 install pyulog``` on MacOS or windows or any other.

(On Mac, it maybe have been installed by the "brew install px4-dev", when you set up a PX4 development environment for macOS. you have to find the path of ulog2csv by "which ulog2csv")

To run:
- Clone or download the repo

- In Matlab run ulog_plot after setup the path and file name

```
WARNING: The scripts ulog2csv, ulog2kml, ulog2rosbag, ulog_extract_gps_dump, ulog_info, ulog_messages, ulog_migratedb and ulog_params are installed in '/Users/mch/Library/Python/3.8/bin' which is not on PATH.
  Consider adding this directory to PATH or, if you prefer to suppress this warning, use --no-warn-script-location.
```
