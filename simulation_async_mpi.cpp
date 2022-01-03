//simulation_async_mpi.cpp
#include <cstdlib>
#include <random>
#include <iostream>
#include <fstream>
#include "contexte.hpp"
#include "individu.hpp"
#include "graphisme/src/SDL2/sdl2.hpp"
#include <chrono>
#include <mpi.h>
void màjStatistique( épidémie::Grille& grille, std::vector<épidémie::Individu> const& individus )//mise à jour statistique
{
    for ( auto& statistique : grille.getStatistiques() )
    {
        statistique.nombre_contaminant_grippé_et_contaminé_par_agent = 0;
        statistique.nombre_contaminant_seulement_contaminé_par_agent = 0;
        statistique.nombre_contaminant_seulement_grippé              = 0;
    }
    auto [largeur,hauteur] = grille.dimension();
    auto& statistiques = grille.getStatistiques();
    for ( auto const& personne : individus )
    {
        auto pos = personne.position();

        std::size_t index = pos.x + pos.y * largeur;
        if (personne.aGrippeContagieuse() )
        {
            if (personne.aAgentPathogèneContagieux())
            {
                statistiques[index].nombre_contaminant_grippé_et_contaminé_par_agent += 1;
            }
            else 
            {
                statistiques[index].nombre_contaminant_seulement_grippé += 1;
            }
        }
        else
        {
            if (personne.aAgentPathogèneContagieux())
            {
                statistiques[index].nombre_contaminant_seulement_contaminé_par_agent += 1;
            }
        }
    }
}
//void afficheSimulation(sdl2::window& écran, épidémie::Grille const& grille, std::size_t jour)
void afficheSimulation(sdl2::window& écran, int largeur_grille,int hauteur_grille,    std::vector<épidémie::Grille::StatistiqueParCase>  statistiques ,int jour)
{

    auto [largeur_écran,hauteur_écran] = écran.dimensions();
    sdl2::font fonte_texte("./graphisme/src/data/Lato-Thin.ttf", 18);
    écran.cls({0x00,0x00,0x00});
    // Affichage de la grille :
    std::uint16_t stepX = largeur_écran/largeur_grille;
    unsigned short stepY = (hauteur_écran-50)/hauteur_grille;
    double factor = 255./15.;

    for ( unsigned short i = 0; i < largeur_grille; ++i )
    {
        for (unsigned short j = 0; j < hauteur_grille; ++j )
        {
            auto const& stat = statistiques[i+j*largeur_grille];
            int valueGrippe = stat.nombre_contaminant_grippé_et_contaminé_par_agent+stat.nombre_contaminant_seulement_grippé;
            int valueAgent  = stat.nombre_contaminant_grippé_et_contaminé_par_agent+stat.nombre_contaminant_seulement_contaminé_par_agent;
            std::uint16_t origx = i*stepX;
            std::uint16_t origy = j*stepY;
            std::uint8_t red = valueGrippe > 0 ? 127+std::uint8_t(std::min(128., 0.5*factor*valueGrippe)) : 0;
            std::uint8_t green = std::uint8_t(std::min(255., factor*valueAgent));
            std::uint8_t blue= std::uint8_t(std::min(255., factor*valueAgent ));
            écran << sdl2::rectangle({origx,origy}, {stepX,stepY}, {red, green,blue}, true);
        }
    }

    écran << sdl2::texte("Carte population grippée", fonte_texte, écran, {0xFF,0xFF,0xFF,0xFF}).at(largeur_écran/2, hauteur_écran-20);
    écran << sdl2::texte(std::string("Jour : ") + std::to_string(jour), fonte_texte, écran, {0xFF,0xFF,0xFF,0xFF}).at(0,hauteur_écran-20);
    écran << sdl2::flush;
}
//void simulation(bool affiche)
void simulation(bool affiche,int rank,int nombre_processeurs, MPI_Comm new_comm)
{
     /**************************************/
     /* create a type for struct Grille */
	    MPI_Datatype mpi_StatistiqueParCase_type;
	    MPI_Type_contiguous(3, MPI_INT, &mpi_StatistiqueParCase_type);
	    MPI_Type_commit(&mpi_StatistiqueParCase_type);
	   
    /**************************************/

    
    constexpr const unsigned int largeur_écran = 1280/2, hauteur_écran = 1024/2;// j'ai changé dimension de l'écran
    sdl2::window écran("Simulation épidémie de grippe", {largeur_écran,hauteur_écran});

    unsigned int graine_aléatoire = 1;
    std::uniform_real_distribution<double> porteur_pathogène(0.,1.);


    épidémie::ContexteGlobal contexte;
    // contexte.déplacement_maximal = 1; <= Si on veut moins de brassage
    contexte.taux_population = 200'000;
    //contexte.taux_population = 1'000;
    contexte.interactions.β = 60.;
    std::vector<épidémie::Individu> population;
    population.reserve(contexte.taux_population);
    épidémie::Grille grille{contexte.taux_population};

    if (rank != 0) // le seul pour l'affichage est le processeur 0 
    {
        	sdl2::finalize();
    }

    auto [largeur_grille,hauteur_grille] = grille.dimension();
    auto statistiques = grille.getStatistiques();
    // L'agent pathogène n'évolue pas et reste donc constant...
    épidémie::AgentPathogène agent(graine_aléatoire++);
    // Initialisation de la population initiale :
    
    for (std::size_t i = 0; i < contexte.taux_population; ++i )
    {
		std::default_random_engine motor(100*(i+1));
		population.emplace_back(graine_aléatoire++, contexte.espérance_de_vie, contexte.déplacement_maximal);
		population.back().setPosition(largeur_grille, hauteur_grille);
		if (porteur_pathogène(motor) < 0.2)
		{
		    population.back().estContaminé(agent);   
		}
    }

    std::size_t jours_écoulés = 0;
    int         jour_apparition_grippe = 0;
    int         nombre_immunisés_grippe= (contexte.taux_population*23)/100;
    sdl2::event_queue queue;

    bool quitting = false;

    std::ofstream output("Courbe.dat");
    output << "# jours_écoulés \t nombreTotalContaminésGrippe \t nombreTotalContaminésAgentPathogène()" << std::endl;

    épidémie::Grippe grippe(0);


    std::cout << "Début boucle épidémie" << std::endl << std::flush;

    while (!quitting)
    {
	    std::cout<<"while here we go"<<std::endl;
	    if (rank>0)
	    {
		std::chrono::time_point < std::chrono::system_clock > start, end;
		start = std::chrono::system_clock::now();
	       //Début temporisation
	       std::cout<<"rank = 1 "<<std::endl;


		auto events = queue.pull_events();
		for ( const auto& e : events)
		{
		    if (e->kind_of_event() == sdl2::event::quit)
		        quitting = true;
		}
		if (jours_écoulés%365 == 0)// Si le premier Octobre (début de l'année pour l'épidémie ;-) )
		{
			    grippe = épidémie::Grippe(jours_écoulés/365);
			    jour_apparition_grippe = grippe.dateCalculImportationGrippe();
			    grippe.calculNouveauTauxTransmission();
			    // 23% des gens sont immunisés. On prend les 23% premiers
			    for ( int ipersonne = 0; ipersonne < nombre_immunisés_grippe; ++ipersonne)
			    {
				population[ipersonne].devientImmuniséGrippe();
			    }
			    for ( int ipersonne = nombre_immunisés_grippe; ipersonne < int(contexte.taux_population); ++ipersonne )
			    {
				population[ipersonne].redevientSensibleGrippe();
			    }
		}
			
		if (jours_écoulés%365 == std::size_t(jour_apparition_grippe))
		{
			    for (int ipersonne = nombre_immunisés_grippe; ipersonne < nombre_immunisés_grippe + 25; ++ipersonne )
			    {
				population[ipersonne].estContaminé(grippe);
			    }
		}
		// Mise à jour des statistiques pour les cases de la grille :
		màjStatistique(grille, population);
			
		// On parcout la population pour voir qui est contaminé et qui ne l'est pas, d'abord pour la grippe puis pour l'agent pathogène
			
		std::size_t compteur_grippe = 0, compteur_agent = 0, mouru = 0;
		
		
		//for ( auto& personne : population )
		for (int i=(rank-1)*population.size()/(nombre_processeurs-1) ; i<rank*population.size()/(nombre_processeurs-1); i++)
													 // ici on fait parallelisme de ce boucle
													 // pour le groupe de processeurs reservés
													 // pour la simulation.
		{
			    if (population[i].testContaminationGrippe(grille, contexte.interactions, grippe, agent))
			    {
				compteur_grippe ++;
				population[i].estContaminé(grippe);
			    }
			    if (population[i].testContaminationAgent(grille, agent))
			    {
				compteur_agent ++;
				population[i].estContaminé(agent);
			    }
			    // On vérifie si il n'y a pas de personne qui veillissent de veillesse et on génère une nouvelle personne si c'est le cas.
			    if (population[i].doitMourir())
			    {
				mouru++;
				unsigned nouvelle_graine = jours_écoulés + population[i].position().x*population[i].position().y;
				population[i] = épidémie::Individu(nouvelle_graine, contexte.espérance_de_vie, contexte.déplacement_maximal);
				population[i].setPosition(largeur_grille, hauteur_grille);
			    }
			    population[i].veillirDUnJour();
			    population[i].seDéplace(grille);
		}
		
		MPI_Gather (population.data(),  
		population.size()*sizeof(épidémie::Individu)/((nombre_processeurs-1)),
		MPI_PACKED, &population[0],  population.size()*sizeof(épidémie::Individu)/((nombre_processeurs-1)), 
		MPI_PACKED, 0, new_comm );
		
		int grippe = compteur_grippe;
            	int agent = compteur_agent;
            	int mor = mouru;
            	
            	MPI_Allreduce(&grippe, &compteur_grippe, 1, MPI_INT, MPI_SUM, new_comm);
            	MPI_Allreduce(&agent, &compteur_agent, 1, MPI_INT, MPI_SUM, new_comm);
            	MPI_Allreduce(&mor, &mouru, 1, MPI_INT, MPI_SUM, new_comm);
		
		
		MPI_Send (grille.getStatistiques().data() , largeur_grille*hauteur_grille, mpi_StatistiqueParCase_type , 0, 0, MPI_COMM_WORLD );
		MPI_Send (&jours_écoulés , 1 , MPI_INT , 0, 0, MPI_COMM_WORLD );
		
		
		// fin de temporisation 
		end = std::chrono::system_clock::now();

		std::chrono::duration < double >elapsed_seconds = end - start;
		std::cout<<"le temps passé dans la simulation par pas de temps : " << elapsed_seconds.count() << " secondes."<<std::endl;
		
               }
		//#########################################################
		//##    Affichage des résultats pour le temps  actuel    ##
		//#########################################################
		else
		if (affiche) 
		{
			int jour;
		    	MPI_Status status ;
			
			MPI_Recv (statistiques.data() , largeur_grille*hauteur_grille , 
			mpi_StatistiqueParCase_type , 1, 0, MPI_COMM_WORLD , & status );
	    		
	    		MPI_Recv (&jour , 1, MPI_INT , 1, 0, MPI_COMM_WORLD ,& status );

			afficheSimulation(écran, largeur_grille,hauteur_grille,statistiques ,jour);
		}
	       /* std::cout << jours_écoulés << "\t" << grille.nombreTotalContaminésGrippe() << "\t"
		          << grille.nombreTotalContaminésAgentPathogène() << std::endl;*/
		output << jours_écoulés << "\t" << grille.nombreTotalContaminésGrippe() << "\t"
		       << grille.nombreTotalContaminésAgentPathogène() << std::endl;
		jours_écoulés += 1;
	}// Fin boucle temporelle
	output.close();
}

int main(int argc, char* argv[])
{
    
	int nombre_processeurs , rank;
	

	MPI_Init ( &argc , &argv );
	
	MPI_Comm_size ( MPI_COMM_WORLD , & nombre_processeurs );

	MPI_Comm_rank ( MPI_COMM_WORLD , &rank);
 
    // Determine the colour and key based on whether my rank is even.
	// ici on travaille avec deux processeur donc on peut affecter 

    int color=rank;
    int key=rank;
    MPI_Comm new_comm;
   
    MPI_Comm_split(MPI_COMM_WORLD, color, key, &new_comm);
    
    
    
    bool affiche = true;
    for (int i=0; i<argc; i++) 
    {
      	std::cout << i << " " << argv[i] << "\n";
    
    	if (std::string("-nw") == argv[i]) affiche = false;
    }
  
    sdl2::init(argc, argv);
    {
        simulation(affiche,rank,nombre_processeurs,new_comm);
    }
    
    sdl2::finalize();
    MPI_Finalize ();
  
  
    return EXIT_SUCCESS;
}
