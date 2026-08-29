This repo contains a pipeline i proposed and created in 10th and 11th grade that generates new candidate molecules for hypertension drugs, and then predicts their binding affinity to AT1R (how good the drug is).
This research was done under the supervision of Dr. Sadashiva Karnik in his lab in the Cleveland Clinic, I wouldn't have got nearly as far as I did without the guidance provided by members of his lab (tysm if any of you are reading this).
All code written in this repository is my own, certain guiding ideas/principles were given to me from lab members.

Ngl there is a pretty decent chunk of stuff missing from here, 
I remember having some trouble pushing stuff to git because of how big the repo got, and since I was the only collaborator I just gave up lol...

Most of the interesting stuff that I wrote and thought of myself is in Binding_Affinity_Prediciton, its random forest regression in order to predict the binding affinity of the proposed ligand. 
The innovative things I did lied within the feature engineering, there was a LOT of it!!

Molecule_Generation houses an implementation of Microsfot MoLeR that I fine-tuned to generate reasonable ligands for the AT1R receptor, trained with ARBs that are currently in circulation (clinically tested, I think I used 8)
