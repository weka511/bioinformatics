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
 
 #include <sstream>
 #include "catch.hpp"
 #include "newick.hpp"
 #include "tokenizer.hpp"
 
 class Displayer: public Parser::Visitor {
	 void accept(Parser::NewickNode* node,const int depth){
		 stringstream leader;
		 for (auto i=0;i<depth;i++)
			 leader << "-";
		 cout << __FILE__ << " " << __LINE__  << ": " <<leader.str() <<*node  << endl;
	 }
 };
 
TEST_CASE( "Newick Tests", "[newick]" ) {
	Parser parser;
	Tokenizer tokenizer;
	Displayer displayer;
	
	SECTION("get_first_comma_at_top_level") {
		auto tokens = tokenizer.tokenize("A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5");	
		Parser::BranchSet branch_set;
		REQUIRE(branch_set.get_first_comma_at_top_level(tokens) == 3);
	}
	
	SECTION("get_first_comma_at_top_level(1)") {
		auto tokens = tokenizer.tokenize("(C:0.3,D:0.4)E:0.5,A:0.1,B:0.2");	
		Parser::BranchSet branch_set;
		REQUIRE(branch_set.get_first_comma_at_top_level(tokens) == 12);
	}

	SECTION("a naked leaf node ") {
		shared_ptr<Parser::Tree> tree = parser.parse("SomeLeaf;");
		shared_ptr<Parser::Visitor> displayer = make_shared<Displayer>();
		tree->descend(displayer);
		shared_ptr<Parser::NewickNode> subtree = tree->get_child(0);
		shared_ptr<Parser::NewickNode> leaf = subtree->get_child(0);
		shared_ptr<Parser::NewickNode> name = leaf->get_child(0);
		REQUIRE(name->get_name() == "SomeLeaf");
	}
	
	SECTION("a naked leaf node with length") {
		auto tokens = tokenizer.tokenize("SomeLeaf:3.1415;");	
		Parser::Tree tree;
		tree.parse(tokens);
		shared_ptr<Parser::NewickNode> subtree = tree.get_child(0);
		shared_ptr<Parser::NewickNode> leaf = subtree->get_child(0);
		shared_ptr<Parser::NewickNode> name = leaf->get_child(0);
		REQUIRE(name->get_name() == "SomeLeaf");
		shared_ptr<Parser::NewickNode> length = name->get_child(0);
		REQUIRE(length->get_length() == 3.1415);
	}
	
	SECTION("Test parser with labels") {
		shared_ptr<Parser::Tree> tree = parser.parse("(A,B,(C,D));");
		shared_ptr<Parser::Visitor> displayer = make_shared<Displayer>();
		tree->descend(displayer);
	}
	
	SECTION("Test parser with labels and lengths(1)") {
		shared_ptr<Parser::Tree> tree = parser.parse("(A:2.1,B:1.2,(C:3.0,D:0.2));");
		shared_ptr<Parser::Visitor> displayer = make_shared<Displayer>();
		tree->descend(displayer);
	}
	
	// SECTION("Test parser with labels") {
		// auto tree = parser.parse(
						// "(((Abantias_enhydris,(Bradypodion_australis,Lasiodora_zenobia))," 
						// "(Bubulcus_carnivorus,Riparia_rostratus)),(((((((Aegialites_chukar,(((Anolis_schokari,(Aphonopelma_galeatus," 
						// "(((Callipogon_fernandi,(Machetes_belliana,Scincus_fuellebornii)),Phelsuma_irregularis),"
						// "(Citellus_carbonaria,Haliaetus_emilia)))),Pareas_acanthinura),((Anthropoides_ferruginea,"
						// "(Grampus_tataricus,Prunella_lineatus)),Burhinus_conicus))),(Hyla_blythi,Tursiops_heliaca)),"
						// "Terpsihone_rufodorsata),((((((Allobates_aegagrus,Perdix_clypeatus),Pachytriton_hilarii),"
						// "Aythia_caudatus),Babycurus_bicinctores),Chrysemys_fasciolata),((Cottus_fernandi,Ctenosaura_vertebralis),"
						// "Minipterus_ciliatus))),Antaresia_bengalensis),Ethmostigmus_pardus),Alaus_tristis),"
						// "(((((((((((((Aegialites_mysticetus,Phyllopneuste_mandarina),Lyrurus_dendrophila),"
						// "((Phyllopneuste_viscivorus,Sphenops_sarasinorum),Xenopeltis_maihensis)),Eublepharis_dignus),"
						// "(Emberiza_dolosus,Petrocincla_dactylisonans)),(Sorex_difficilis,Spermophilus_stejnegeri)),"
						// "(Dipsosaurus_meermani,Testudo_lividum)),((Furcifer_horridum,Tropidurus_wogura),"
						// "(Middendorffinaia_maihensis,Motacilla_nigropalmatus))),Minipterus_cranwelli),"
						// "((((((((((((((Ahaetulla_tetrax,((((Carabus_americanus,(Haplopelma_bicoloratum,Leuciscus_himalayanus)),"
						// "Oxyura_sirtalis),((Chlamydotis_hendersoni,Lanius_gibbosus),Siniperca_lividum)),Coenobita_borealis)),"
						// "(Archispirostreptus_linaria,Phalacrocorax_vitulina)),Apodora_cinereus),Pseudorca_geniculata),Aythya_aleutica),"
						// "((((((((Allobates_turtor,Vulpes_taeniura),Psammophis_alcinous),Himantopus_orientalis),((Anthropoidae_kurilensis,"
						// "Pyxicephalus_cenchria),Bombyx_saxatilis)),Rhamphiophis_ferox),(((Apus_fluviatilis,Brachypelma_marmorata),"
						// "Ctenotus_means),Otis_gregarius)),Elaphe_coturnix),Calotes_milii)),(Bombus_glacialis,"
						// "(Dryobates_rutila,Holaspis_pachypus))),Bubulcus_apus),Lampropeltis_proteus),Pterocles_gallinago),"
						// "Phoca_naumanni),(Eublepharis_leucoptera,Pyrgilauda_corone)),Litoria_parahybana),"
						// "((Budytes_rutila,Sceloporus_constricticollis),Phylloscopus_pallasii))),Tryngites_cliffordii),"
						// "Atrophaneura_lagopus),Rhinolophus_sanguinolentus));");
		// auto displayer = make_shared<Displayer>();
		// tree->descend(displayer);
	// }
	
 }