# Aligner
This project was my research work in the company EPAM Systems. 
This script searches for the optimal offset for the pairwise alignment of two restriction maps for the given parameters of the error (skip the cut and the error in the length of the fragment).

It was implemented two approaches:
* Offset in fragments
![alt text](https://pp.vk.me/c630529/v630529883/f9e8/GYggrlzge0w.jpg "Fragment approach")
* Offset in cuts
![alt text](https://pp.vk.me/c630529/v630529883/f9f1/GHQfEzuYEmA.jpg "Kilobase approach")

# Algorithms

* **Offset in fragments**:
Since the offset is an integer, then iterate over all of the offset generated and for each pair of restriction maps are building the region and consider how many times we were wrong

* **Offset in kilobase**:
Since the offset is an double, we use binary search to find optimal offset.

# Results
![alt text](https://pp.vk.me/c630529/v630529883/f9d4/MkpicClwjTQ.jpg "Fragment approach")
![alt text](https://pp.vk.me/c630529/v630529883/f9de/jxKP18Wmalo.jpg "Kilobase approach")
