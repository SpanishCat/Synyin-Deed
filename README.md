# Synyin Deed

# General
Synteny Index (SI) is a newly-developed method for genetic alignment of organisms. It has been proved to be more accurate (in certain conditions) than sequence-based approaches. That is, for prokaryotes (Bacteria & Archea).
I've been developing Syntin Deed as part of a study, which undermines to find out whether SI is also as good for eukaryotes (All other organisms) genetic comparisons.

# How to use?
Download the program to your desktop and run "run.bat".
You should see something like this:
![image](https://user-images.githubusercontent.com/73846269/204306416-6f4aa33d-788f-4668-8efd-a5bb5990e801.png)


"Covnert" and "Analyze" are essentials. The rest are quality-of-life functions.
Press the number of category of commands that you want to choose. For example: Enter '2' to enter the "Convert" category.


1. Multiples:
Allows you to convert/inspect all of the files in input folder one after the other.

2. Convert:

![image](https://user-images.githubusercontent.com/73846269/204307786-125b1b04-93a2-47c3-bad6-d4785434b947.png)

1. GFF3 & FASTA -> GENELIST2: This is the first essential command. A .genelist2 file is a special type of file that you can then apply the SI algorithm to.
This function takes 4 files: The GFF and FASTA of each genome. Put those of the first genome in Main_Input/File 1 and those of the second one in Main_Input/File 2, then choose this command in console.

Output:
All output files are generated in Main_Output folder.


# Developer Warning
###Read the following if you intend on inspecting/modifying the source code:
Because of the typical sheer size of genomes entered into the program, which it needs to go over numerous times, every little inefficiency in many areas of the code can easily add a few more minutes to runtime. 
Over the course of the entire code this can result in runtime many times the current one. 

Because of that, this code makes an extensive use of tuples, tuple comprehensions (tuple(x for x in < iterbale >)) and multiprocess-mapping.
This sometimes reduce readability, especially when tuple comprehension is used instead of for loop, making them a lot faster but also harder to read.

**You were warned!**

# Further Reading
https://www.sciencedirect.com/science/article/pii/S1055790318303087?casa_token=by3jqQtFBGkAAAAA:3cTVYapeSlpMskMnUAtWrLXgZej8BnTzTVLnqRDD8rQLI_OyPohnS4gVz2NGjmJP5evNGClh3A
