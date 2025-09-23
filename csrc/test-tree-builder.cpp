/**
 * Copyright (C) 2025 Simon Crase
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this software.  If not, see <http://www.gnu.org/licenses/>
 *
 * (,,(,));                               no nodes are named
 * (A,B,(C,D));                           leaf nodes are named
 * (A,B,(C,D)E)F;                         all nodes are named
 * (:0.1,:0.2,(:0.3,:0.4):0.5);           all but root node have a distance to parent
 * (:0.1,:0.2,(:0.3,:0.4):0.5):0.0;       all have a distance to parent
 * (A:0.1,B:0.2,(C:0.3,D:0.4):0.5);       distances and leaf names (popular)
 * (A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;     distances and all names
 * ((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;    a tree rooted on a leaf node (rare)
 */
 
 #include "catch.hpp"

 #include "tree-builder.hpp"
 
 class Displayer: public Node::Visitor {
	 void accept(Node& node,const int depth){
		 stringstream leader;
		 for (auto i=0;i<depth;i++)
			 leader << "-";
		 cout << __FILE__ << " " << __LINE__  << ": " /*<<leader.str()*/ <<node  << endl;
	 }
 };
 

 
 TEST_CASE( "Tree Builder Tests", "[tree_builder]" ) {
	Parser parser;
	shared_ptr<TreeBuilder> builder = make_shared<TreeBuilder>();
	Displayer displayer;
	Node::reset();
	
	SECTION("Test parser with labels") {
			shared_ptr<Parser::Tree> tree = parser.parse("(A,B,(C,D));");
			tree->descend(builder);
			shared_ptr<Node> new_tree = builder->get_result();
			new_tree->visit(displayer);
	}
	 
	SECTION("Test parser with labels and lengths(1)") {
		auto new_tree = TreeBuilder::create("(A:2.1,B:1.2,(C:3.0,D:0.2));");
		new_tree->visit(displayer);
	}
	
	SECTION("Test parser with labels") {
		auto tree = TreeBuilder::create(
						"(((Abantias_enhydris,(Bradypodion_australis,Lasiodora_zenobia))," 
						"(Bubulcus_carnivorus,Riparia_rostratus)),(((((((Aegialites_chukar,(((Anolis_schokari,(Aphonopelma_galeatus," 
						"(((Callipogon_fernandi,(Machetes_belliana,Scincus_fuellebornii)),Phelsuma_irregularis),"
						"(Citellus_carbonaria,Haliaetus_emilia)))),Pareas_acanthinura),((Anthropoides_ferruginea,"
						"(Grampus_tataricus,Prunella_lineatus)),Burhinus_conicus))),(Hyla_blythi,Tursiops_heliaca)),"
						"Terpsihone_rufodorsata),((((((Allobates_aegagrus,Perdix_clypeatus),Pachytriton_hilarii),"
						"Aythia_caudatus),Babycurus_bicinctores),Chrysemys_fasciolata),((Cottus_fernandi,Ctenosaura_vertebralis),"
						"Minipterus_ciliatus))),Antaresia_bengalensis),Ethmostigmus_pardus),Alaus_tristis),"
						"(((((((((((((Aegialites_mysticetus,Phyllopneuste_mandarina),Lyrurus_dendrophila),"
						"((Phyllopneuste_viscivorus,Sphenops_sarasinorum),Xenopeltis_maihensis)),Eublepharis_dignus),"
						"(Emberiza_dolosus,Petrocincla_dactylisonans)),(Sorex_difficilis,Spermophilus_stejnegeri)),"
						"(Dipsosaurus_meermani,Testudo_lividum)),((Furcifer_horridum,Tropidurus_wogura),"
						"(Middendorffinaia_maihensis,Motacilla_nigropalmatus))),Minipterus_cranwelli),"
						"((((((((((((((Ahaetulla_tetrax,((((Carabus_americanus,(Haplopelma_bicoloratum,Leuciscus_himalayanus)),"
						"Oxyura_sirtalis),((Chlamydotis_hendersoni,Lanius_gibbosus),Siniperca_lividum)),Coenobita_borealis)),"
						"(Archispirostreptus_linaria,Phalacrocorax_vitulina)),Apodora_cinereus),Pseudorca_geniculata),Aythya_aleutica),"
						"((((((((Allobates_turtor,Vulpes_taeniura),Psammophis_alcinous),Himantopus_orientalis),((Anthropoidae_kurilensis,"
						"Pyxicephalus_cenchria),Bombyx_saxatilis)),Rhamphiophis_ferox),(((Apus_fluviatilis,Brachypelma_marmorata),"
						"Ctenotus_means),Otis_gregarius)),Elaphe_coturnix),Calotes_milii)),(Bombus_glacialis,"
						"(Dryobates_rutila,Holaspis_pachypus))),Bubulcus_apus),Lampropeltis_proteus),Pterocles_gallinago),"
						"Phoca_naumanni),(Eublepharis_leucoptera,Pyrgilauda_corone)),Litoria_parahybana),"
						"((Budytes_rutila,Sceloporus_constricticollis),Phylloscopus_pallasii))),Tryngites_cliffordii),"
						"Atrophaneura_lagopus),Rhinolophus_sanguinolentus));");
		tree->visit(displayer);
	}
}