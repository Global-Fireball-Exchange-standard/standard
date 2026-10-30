# The Global Fireball Exchange (GFE) standard

This repository contains documentation of the standard for exchange of astrometric and photometric data between fireball or meteor camera systems.  The data relates to the capture by a single meteor camera of a single meteor path.  It is intended to include all trajectory data in a single data file, including for multiple fragments and for different positional measurements. 

Sample data files in GFE format are included in the repository, and correspond to actual observations from five different camera systems of the fall of the Winchcombe meteorite in the UK on 2021-02-28.  In all, sixteen observations were made from five different camera systems, though only five sample files have been uploaded here - one from each type of system.

Jupyter code to read files in GFE format and to generate a sample GFE data file is also included, as is a Jupyter script which converts fireball single-station observations between formats.

To check whether a file conforms to the standard, use the online GFE format linter at https://global-fireball-exchange-standard.github.io/standard/ (it runs in your browser; files are not uploaded), or run `GFE_format_linter.ipynb` in Jupyter. Both use the checks in `gfe_linter.py`.

How did we come up with the GFE format?  It is very closely based on the Desert Fireball Network's data file format.  The selection process was completed in 2020 and is detailed in a longform draft paper available in this repository, as well as in a shortform slide pack, and also in a 9.5 minute video which is available here: https://player.vimeo.com/video/458077336
