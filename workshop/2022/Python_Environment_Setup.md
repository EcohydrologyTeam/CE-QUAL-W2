# CE-QUAL-W2 Workshop: Python and JupyterLab

The following instructions will guide you in how to use the Anaconda distribution of Python to:

1. Install a custom library of software you will need for this workshop
2. Run JupyterLab, a modern computing environment that allows you to run Python code interactively. It also contains a number of useful features, such as a file browser and image and PDF viewer.

Once you have JupyterLab running, you will be able to open and run one of the provided "Jupyter Notebook" files.

##### 1. Open your Anaconda Prompt. This is also known as a terminal. It will allow you to enter a series of commands. These will be shown below.

##### 2. Change to the directory containing the files we will be working with for workshop 2.11. This will ensure that everything you need for the workshop will be easy to access. You will first need to change directories within the Anaconda prompt. Please do the following:

* Type the letters "cd" without the quotes. This stands for "change directory".
* Type one space on the spacebar.
* Go to your file browser (Windows Explorer). Drag and drop the folder that contains the 2.11_Workshop_Introduction_to_Python.ipynb file to the Anaconda Prompt.
* Verify that you have changed to the right directory.

##### 3. Create a new virtual environment, a sandbox that will contain only the software you want to use. Please type or copy-paste the following in the Anaconda Prompt:

conda create --name w2 python=3.9

##### 4. Activate your new environment. Enter the following:

conda activate w2

##### 5. Install the software libraries will will require for this workshop. Enter the following:

conda install numpy matplotlib seaborn pandas h5py yaml bokeh ipython openpyxl jupyter jupyterlab

##### 6. If the above steps have run correctly, without errors. You are now ready to run JupyterLab. Please enter the following. Note: "jupyter lab" is written as two words.

jupyter lab

When JupyterLab launches, you are ready to open a notebook.