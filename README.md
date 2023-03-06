This is a tiny plot tool of matlab for ulog form PX4, which inspired by [matulog](https://github.com/CarlOlsson/matulog). You need to plot what you want by your shelf.

it have been test on MacOS Monterey 12.3, matlab 2021b and PX4 v1.12.3 version.

Requirements: matlab2016a (or a newer version) and pyulog is needed for `
ulog_matlab_plot` to run, install the latter using ```pip install pyulog``` or ```pip3 install pyulog``` on MacOS(it maybe have been installed by the "brew install px4-dev", when you set up a PX4 development environment for macOS. you have to find the path of ulog2csv by "which ulog2csv")

To run:
- Clone or download the repo

- In Matlab run ulog_plot after setup the path and file name
