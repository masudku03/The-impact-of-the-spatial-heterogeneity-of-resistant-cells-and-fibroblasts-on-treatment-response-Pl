# The impact of the spatial heterogeneity of resistant cells and fibroblasts on treatment response
Abstract

A long-standing practice in the treatment of cancer is that of hitting hard with the maximum tolerated dose to eradicate tumors. This continuous therapy, however, selects for resistant cells, leading to the failure of the treatment. A different type of treatment strategy, adaptive therapy, has recently been shown to have a degree of success in both preclinical xenograft experiments and clinical trials. Adaptive therapy is used to maintain a tumor's volume by exploiting the competition between drug-sensitive and drug-resistant cells with minimum effective drug doses or timed drug holidays. To further understand the role of competition in the outcomes of adaptive therapy, we developed a 2D on-lattice agent-based model. Our simulations show that the superiority of the adaptive strategy over continuous therapy depends on the local competition shaped by the spatial distribution of resistant cells. Intratumor competition can also be affected by fibroblasts, which produce microenvironmental factors that promote cancer cell growth. To this end, we simulated the impact of different fibroblast distributions on treatment outcomes.  As proof of principle, we focused on five types of distribution of fibroblasts characterized by different locations, shapes, and orientations of the fibroblast region with respect to the resistant cells. Our simulation shows that the spatial architecture of fibroblasts modulates tumor progression in both continuous and adaptive therapy. Finally, as a proof of concept, we simulated the outcomes of adaptive therapy of a virtual patient with four metastatic sites composed of different spatial distributions of fibroblasts and drug-resistant cell populations. Our simulation highlights the importance of undetected metastatic lesions on adaptive therapy outcomes. 

Articel preprint: https://www.biorxiv.org/content/10.1101/2021.06.01.446525v1

The code ("ProstateCancerGrid.java") is developed using HAL platform(https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007635). It can be used to reproduce the single metastasis results in the above prerpint. To use the code the following steps would work.

1. Install HAL(https://github.com/MathOnco/HAL).
2. Make a folder "ProstateCancer" in the directory "C:\HAL-master\".
3. Copy all the files in the directory "C:\HAL-master\ProstateCancer".
4. Provide one of the csv files as input file in line 192.
5. Set the output file name and directory in line 223.
6. Set the parameter values in lines 123-133 and 180-181.
7. Choose any of the fibroblast strucutres from lines 149-156.
8. Run "ProstateCancerGrid.java"

For any further query please feel free to contact us.
