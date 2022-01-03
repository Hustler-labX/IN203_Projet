# IN203_Projet

J'ai essayé de transformer le projet (en code séquentiel) en plusieurs versions en code parallèle. Donc on a ici 6 versions :

        • simulation_sync_affiche_mpi.cpp : juste une simple implémentation avec MPI_Send / MPI_Recv.
        • simulation_async_affiche_mpi.cpp : à l'aide de MPI_IProbe / MPI_Recv / MPI_ Isend.
        • simulation_async_omp.cpp : utilisation de bibliothèque OpenMP en parallélisant les boucles FOR : 
        en fonction " màjStatistique " et en partie simulation avec la commande " # pragma omp parallel for reduction(+:compteur_grippe,compteur_agent,mouru) ".

        • simulation_omp(personel).cpp : j'ai pensé à une autre solution et je ne suis pas sure si c'est vraiment correct.
        • simulation_async_mpi.cpp : Utilisation de MPI_Split / MPI_Gather / MPI_Allreduce.
        • simulation_async_mpi_omp.cpp : combiner les deux dernières versions.



À chaque compilation, on doit faire des changements en ficher Makefile : soit on change le nom de ficher à compiler (selon la version qu'on souhait avoir) et les include de ficher (Si on travaille avec OpenMP, ça sera avec Make_linux.inc sinon (avec MPI) on a le ficher Make_linux2.inc.

Voir ficher " rapport.opt " pour avoir les résultats.





