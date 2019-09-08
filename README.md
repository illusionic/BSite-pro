# Paper abstract / Intro:

Accurate protein binding site annotations are vital for a profound understanding of biological processes and
protein interactions. Due to less supporting information a large number of proteins remain uncharacterised.
For large sets of uncharacterised proteins, only amino acid information is available. In this paper, we proposed
BSite-pro – a traditional approach – which utilizes only the protein sequences to predict its associated binding
site. The prediction process requires hand-crafted features of the protein sequences. Our results indicate
significant improvements in terms of prediction accuracy and recall when compared to other sequence-based
methods. BSite-pro achieves an overall validation accuracy of 85.06% and recall of 82.17%. Finally, we
discuss that using the same feature extraction methods and model, we successfully solved two different types
of problem i.e. protein active and conserved sites prediction.

# Authors:

* Yasrub Bashir Malik, Khalid Muneer (p146116@nu.edu.pk, p136145@nu.edu.pk) -- Queries about ML should go here.
* Hafeez ur Rehman (hafeez.urrehman@nu.edu.pk) -- Queries about Bioinformatics should go here.

Preprint of related publication available here: https://www.researchgate.net/publication/330556523_BSite-pro_Traditional_Approach_for_Binding_Site_Prediction_in_Protein_Sequence

![Paper Image 1](https://github.com/illusionic/BSite-pro/blob/master/src/interpro-annotations.PNG)

# Import points:
- Requires python3.6
- See requirements.txt for exact version of libraries used.
- It's suggested that you use `virtualenv` to create a new environment and then install required packages(Updated).

    ```
    pip3 install virtualenv
    virtualenv bs
    cd bs
    . bin/activate
    git cone <git_repo_url>
    pip install -r <git_repo_name>src/requirements.txt
    ``` 
    
# Execution:
1) The code is properly commented with guidelines on how to extract features of protein sequence data. 
2) Once the files have been created, one can run train/validation:

This is done through `python train.py`.

# License

This code is provided under the MIT License.

Copyright 2019 Hafeez Ur Rehman and Yasrub Bashir Malik

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
