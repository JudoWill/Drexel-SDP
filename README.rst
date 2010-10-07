This is the beginnings of a GitHub repository for a Bioinformatics based Senior Design Project.


CODING TODO
===========
In no particular order:

- Reformat the code so its a set of functions instead of scripts
  - Probably should get rid of the "folder" idea since it makes it difficult to track information through the function.
  - Probably should create a "main" function which calls the subfunctions and handles the data moving back and forth.
- Clean up code!
  - Remove all EVAL statements.
  - Remove all CLEAR statements.
  - Change all LOAD functions to load into a `struct` instead of loading into the local scope ... this helps make it easy to see where information is comming from.
  - Change any path-related function to use `fullfile` instead of manually constructing the file-path.


CODING EXPERIMENTS TODO
=======================
In no particular order

- Can you find the ruler using only image segmentation?
  - Are they always the same color?
- Can you use OCR to read the 'inches'/'cm' and/or the date?
  - This is probably contingent on the ability to isolate the ruler automatically
- Can you use clustering methods to find the likely color of the wound ... at least from the image I have now it looks like the wound is the most consistent color on the foot.
- Is there a way to automatically isolate a Region of Interest (ROI)? ie. crop out the foot from the exam table, and other background.
- Setup a "report" function using the Matlab publishing tool.
  - You'll need to ask questions about what the users want to actually see in a report.

