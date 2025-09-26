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
#include "newick.hpp"
#include "tokenizer.hpp"
 
  
TEST_CASE( "Newick Tests", "[newick]" ) {
	Parser parser;
	Tokenizer tokenizer;
	
	SECTION("get_first_comma_at_top_level") {
		auto tokens = tokenizer.tokenize("A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5;");	
		Parser::BranchSet branch_set;
		REQUIRE(branch_set.get_first_comma_at_top_level(tokens) == 3);
	}
	
	SECTION("get_first_comma_at_top_level(1)") {
		auto tokens = tokenizer.tokenize("(C:0.3,D:0.4)E:0.5,A:0.1,B:0.2;");	
		Parser::BranchSet branch_set;
		REQUIRE(branch_set.get_first_comma_at_top_level(tokens) == 12);
	}

	SECTION("a naked leaf node ") {
		shared_ptr<Parser::Tree> tree = parser.parse("SomeLeaf;");
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
		auto tree = parser.parse("(A,B,(C,D));");
		auto subtree = tree->get_child(0);
		REQUIRE(subtree->get_type() == Parser::Type::SubTree);
		auto internal = subtree->get_child(0);
		REQUIRE(internal->get_type() == Parser::Type::Internal);
		auto branchset1 = internal->get_child(0);
		REQUIRE(branchset1->get_type() == Parser::Type::BranchSet);
		auto branch1 = branchset1->get_child(0);
		REQUIRE(branch1->get_type() == Parser::Type::Branch);
		auto subtree1 = branch1->get_child(0);
		REQUIRE(subtree1->get_type() == Parser::Type::SubTree);	
		auto leaf1 = subtree1->get_child(0);
		REQUIRE(leaf1->get_type() == Parser::Type::Leaf);
		auto name1 = leaf1->get_child(0);
		REQUIRE(name1->get_type() == Parser::Type::Name);
		REQUIRE(name1->get_name() == "A");
		auto branchset2 = branchset1->get_child(1);
		REQUIRE(branchset2->get_type() == Parser::Type::BranchSet);
		auto branch2 = branchset2->get_child(0);
		REQUIRE(branch2->get_type() == Parser::Type::Branch);
		auto subtree2 = branch2->get_child(0);
		REQUIRE(subtree2->get_type() == Parser::Type::SubTree);	
		auto leaf2 = subtree2->get_child(0);
		REQUIRE(leaf2->get_type() == Parser::Type::Leaf);
		auto name2 = leaf2->get_child(0);
		REQUIRE(name2->get_type() == Parser::Type::Name);
		REQUIRE(name2->get_name() == "B");
		auto branchset3 = branchset2->get_child(1);
		REQUIRE(branchset3->get_type() == Parser::Type::BranchSet);
		auto branch3 = branchset3->get_child(0);
		REQUIRE(branch3->get_type() == Parser::Type::Branch);
		auto subtree3 = branch3->get_child(0);
		REQUIRE(subtree3->get_type() == Parser::Type::SubTree);	
		auto internal2 = subtree3->get_child(0);
		REQUIRE(internal2->get_type() == Parser::Type::Internal);
		auto branchset4 = internal2->get_child(0);
		REQUIRE(branchset4->get_type() == Parser::Type::BranchSet);
		auto branch4 = branchset4->get_child(0);
		REQUIRE(branch4->get_type() == Parser::Type::Branch);
		auto subtree4 = branch4->get_child(0);
		REQUIRE(subtree4->get_type() == Parser::Type::SubTree);	
		auto leaf3 = subtree4->get_child(0);
		REQUIRE(leaf3->get_type() == Parser::Type::Leaf);
		auto name3 = leaf3->get_child(0);
		REQUIRE(name3->get_type() == Parser::Type::Name);
		REQUIRE(name3->get_name() == "C");
		auto branchset5 = branchset4->get_child(1);
		REQUIRE(branchset5->get_type() == Parser::Type::BranchSet);
		auto branch5 = branchset5->get_child(0);
		REQUIRE(branch5->get_type() == Parser::Type::Branch);
		auto subtree5 = branch5->get_child(0);
		REQUIRE(subtree5->get_type() == Parser::Type::SubTree);	
		auto leaf4 = subtree5->get_child(0);
		REQUIRE(leaf4->get_type() == Parser::Type::Leaf);
		auto name4 = leaf4->get_child(0);
		REQUIRE(name4->get_type() == Parser::Type::Name);
		REQUIRE(name4->get_name() == "D");
	}
	
	SECTION("Test parser with labels and lengths(1)") {
		auto tree = parser.parse("(A:2.1,B:1.2,(C:3.0,D:0.2));");
		auto subtree = tree->get_child(0);
		auto internal = subtree->get_child(0);
		auto branchset1 = internal->get_child(0);
		auto branch1 = branchset1->get_child(0);
		auto subtree1 = branch1->get_child(0);
		auto leaf1 = subtree1->get_child(0);
		auto name1 = leaf1->get_child(0);
		REQUIRE(name1->get_name() == "A");
		auto distance1 = name1->get_child(0);
		REQUIRE(distance1->get_length() == 2.1);
		auto branchset2 = branchset1->get_child(1);
		auto branch2 = branchset2->get_child(0);
		auto subtree2 = branch2->get_child(0);
		auto leaf2 = subtree2->get_child(0);
		auto name2 = leaf2->get_child(0);
		REQUIRE(name2->get_name() == "B");
		auto distance2 = name2->get_child(0);
		REQUIRE(distance2->get_length() == 1.2);
		auto branchset3 = branchset2->get_child(1);
		auto branch3 = branchset3->get_child(0);
		auto subtree3 = branch3->get_child(0);
		auto internal2 = subtree3->get_child(0);
		auto branchset4 = internal2->get_child(0);
		auto branch4 = branchset4->get_child(0);
		auto subtree4 = branch4->get_child(0);
		auto leaf3 = subtree4->get_child(0);
		auto name3 = leaf3->get_child(0);
		REQUIRE(name3->get_name() == "C");
		auto distance3 = name3->get_child(0);
		REQUIRE(distance3->get_length() == 3.0);
		auto branchset5 = branchset4->get_child(1);
		auto branch5 = branchset5->get_child(0);
		auto subtree5 = branch5->get_child(0);
		auto leaf4 = subtree5->get_child(0);
		auto name4 = leaf4->get_child(0);
		REQUIRE(name4->get_name() == "D");
		auto distance4 = name4->get_child(0);
		REQUIRE(distance4->get_length() == 0.2);
	}
	
	SECTION("Test parser with QRTD data") {
		auto tree = parser.parse(
								"((Abantias_cioides,((Chlamydotis_hypoleucos,Onychodactylus_kraepelini),Phelsuma_pulchripes)),"
								"((((((((((((((((((((((((((((((((((((((((((((((((((((((((Abantias_horridum,(((((((((((((Agama_sanguinolentus,"
								"Limosa_limosa),Oenanthe_fasciolata),Chrttusia_tinctorius),Mylopharyngodon_leucorodia)"
								",((Ketupa_pedo,Odonthurus_pendulinus),Thecla_trianguligerus)),Otocoris_guineti),((((Alcedo_rubida,"
								"Psammophis_maritimus),Laudakia_adamsii),Ursus_cinerea),((Aphonopelma_avinivi,(Atrophaneura_aleutica"
								",Micropalama_hypomelus)),Buthacus_hyperboreus))),(Bufo_ochropus,Opheodrys_perdix)),Melanocorypha_nigrolineatus),"
								"Pterinochilus_armata),Megaloperdix_javanica),(((Alloporus_cyanea,Telescopus_ameiva),Eremophila_orientalis),"
								"Middendorffinaia_elaphus)),(Spalax_collaris,Totanus_scripta))),Norops_cygnoides),(Amphiuma_barbata,"
								"Gongylophis_ornata)),(((Alectoris_tataricus,Apus_longipes),Himantopus_caudatus),(((((((((Antaresia_circia,"
								"((((((Aphonopelma_multituberculatus,((((Babycurus_falcinellus,Physignathus_laticauda),Clemmys_subrubrum),"
								"Erpeton_sujfunensis),(Bombyx_varius,Erpeton_lineatus))),Eudrornias_himalayanus),(Butastur_hypomelus,"
								"Porzana_tuberculosus)),Fregilegus_quinquestriatus),Cardiocranius_hardwickii),((Circus_aristotelis,"
								"Motacilla_radiata),(Eunectes_calidris,Scolopendra_insularis)))),(Buthacus_bukhunensis,Cyriopagopus_minor))"
								",Nerodia_deremensis),Atrophaneura_ferina),(Paradoxornis_mackloti,Phrynosoma_armata)),(((((((((Asthenodipsas_lehmanni"
								",(Eublepharis_diadema,Pagophila_franckii)),Melanocorypha_guangxiensis),Totanus_exquisita),Lyrurus_clarus),"
								"Sorex_pedo),(Brachyramphus_oenanthe,(Ceratophrys_pholeter,Kinixys_leucoryphus))),Ninox_flammea),"
								"(Rhynchaspis_marinus,Tetraogallus_onocrotalus)),Vormela_solitaria)),Elaphe_plumifrons),Hyperoodon_leucophyllata),"
								"(((Bombina_Jankowskii,Rhynchophis_monoceros),Lyrurus_molurus),Lasiodora_cliffordii)))),Elaphe_fernandi),"
								"(((((((((((Anthropoidae_smaragdina,Leptopelis_baibacina),Boiga_avicularia),Phasianus_leucostomum),"
								"((Canis_irregularis,Platalea_ridibundus),Uromastyx_cinerea)),Ciconia_squaterola),(Eudramias_mykiss,"
								"Furcifer_ceterus)),(Lanius_boulengeri,Monachus_arenarius)),Holaspis_leucogaster),(Notophthalmus_jaculus,"
								"(Rhombomys_armeniacus,Tupinambus_occitanus))),Erpeton_temminckii),(Pterinochilus_murinus,Rhabdophis_hosii))),"
								"Odonthurus_bukhunensis),(((((((((((((((((((Acanthis_avinivi,((Caiman_rosea,(Chelydra_punctatus,Plethodon_limosa))"
								",Phalacrocorax_longipennis)),Onychodactylus_metallica),((Acanthis_plathyrhychos,Pogona_rosmarus),Spalax_czerskii)),"
								"Tiliqua_pulchra),Mogera_miliaris),Porzana_moschata),Almo_falcipennis),Leiurus_not),Natriciteres_monachus),"
								"Ortigometra_ichthyaetus),((((((((((((((((((((((((((Acanthogonatus_mongolica,((Monachus_dactylisonans,"
								"Thecla_zenobia),Salmo_cristatella)),((Gambelia_livia,Psalmopoeus_gebleri),Melanocorypha_glacialis))"
								",Latastia_nasicus),Motacilla_citrsola),Gypaetus_leporosum),(Pachydactylus_naumanni,Psalmopoeus_guineti))"
								",Gyps_vitticeps),Kassina_collybitus),(Butastur_arcticus,((((Ceratophrys_grandis,Gecarcinus_temminckii),"
								"Enhydra_calvus),Nhandu_auriculatus),Homalopsis_pulchra))),(Cygnus_pictus,Gerrhosaurus_dominicus)),"
								"Circaetus_smaragdina),Trapelus_hongkongensis),((Alpes_cristatella,(Enhudra_kopsteini,Liasis_obsoleta)),"
								"Tryngites_torquata)),Pterocles_caudicinctus),Neophron_albopillosum),(Machetes_hasselquistii,Podoces_argentatus))"
								",Citellus_molurus),((((Anas_difficilis,Pelomedusa_caudatus),(Lamprolepis_casualis,Pedostibes_imperator)),"
								"Prunella_sagrei),Ursus_turtor)),(Aegialifes_australis,Salamandra_leucopsis)),((Aegialifes_pulchripes,"
								"((Dryobates_melanoleucus,Pareas_amethistina),Gecarcinus_unicolor)),Callipogon_leiosoma)),(Cuora_albatrus,"
								"Enhudra_pugnax)),Euspiza_flavirufa),((((((((Acheron_mongolica,Eirenis_taezanowskyi),Riparia_fuliginosus),"
								"Streptopelia_graculus),Alectoris_argentatus),(Eudramias_eximia,Squaterola_blythi)),Ingerophrynus_torquatus),"
								"Upupa_monoceros),(Iomachus_clericalis,Uromastyx_tataricus))),Odonthurus_ovata),((((Dendrelaphis_tinctorius,"
								"Salvelinus_tadorna),Epicrates_conicus),Thecla_parahybana),((Lepus_ibera,Myotis_fulva),Trapelus_stagnalis))),"
								"Ardea_deremensis)),Osmoderma_pallidus),Triturus_querquedula),Cynops_cambridgei),Sitta_gallinago),"
								"((Eulabeia_alpestris,Seokia_auriculatus),Megaptera_flavescens)),Vanellus_avocetta),Xenopeltis_dentatus),"
								"((((Aegialites_exanthematicus,Turdus_chrysargos),((((((((Alloporus_oenanthe,(((((((((Chondropython_soloensis,"
								"(Hadogenes_fernandi,Trapelus_constricticollis)),(Grus_metallica,Procellaria_vitulina)),(((Lamprophis_argentatus,"
								"Madagascarophis_dennysii),Sus_semipalmatus),Salvelinus_dubius)),Gekko_rapax),Cinclus_cachinans),"
								"Perdix_tatarica),(((Coleonyx_picta,Podoces_brandtii),Phormictopus_sirtalis),Perdix_tarandus)),"
								"Mogera_ussuriensis),Cuculus_rosmarus)),(Cardiocranius_ceterus,Ptyodactylus_rutilans)),"
								"(Scaphiophryne_not,Scaphiopus_leucomelas)),((((((Castor_carbonaria,((Furcifer_avicularia,Tadarida_cenchria),"
								"Natriciteres_truncatus)),Mabuya_clypeata),Colaeus_caudicinctus),Chalcides_sieboldii),(Kinixys_deremensis,"
								"Regulus_sarasinorum)),(Remiz_kopsteini,Seokia_madagascariensis))),((Chelus_avinivi,Grammostola_pygmeus),"
								"Gyps_ochropus)),Ciconia_orientalis),Dafila_quadrivirgata),(Anodonta_buccata,Capella_rubida))),Litoria_adspersus),"
								"((((Capra_dexter,((((Ciconia_amboinensis,Latastia_davidiana),Lampropeltis_collaris),Homopholis_leucopsis),"
								"(Mesoplodon_caudatus,Phrynohyas_ussuriensis))),Chrttusia_merganser),Rhacophorus_gebleri),"
								"Eucratoscelus_tentaculatum)))),Rhynchophis_dactylisonans),Chelodina_rubicola),(Alopex_leucocephala,"
								"(Butastur_aspera,Pelomedusa_cyanus))),(Capra_pachypus,Elseya_caniceps)),(Buthacus_vipio,"
								"Heterodon_ignicapillus)),Epicrates_schokari),(((((((((((((((((((((Acanthoceros_ampullatus,Oxyura_diffidens),"
								"((((((((((Apodora_grunniens,(Saxicola_aceras,Tringa_alpestris)),Triturus_multituberculatus),"
								"(Litoria_insignis,Teratoscincus_hypoleucus)),Damon_sebae),(((((Boa_enhydris,Platemys_Bernicla),"
								"(((Chondropython_capra,Circaetus_scutulata),Pyrgilauda_unicolor),Procellaria_coelestis)),"
								"Podoces_caesius),Holodactylus_helena),Glareola_torquatus)),(Hottentotta_maritimus,crecca_gregaria)),"
								"Elseya_chrysaetos),Margaritifera_temminckii),Citellus_calidris),Rhabdophis_pallidus)),((Aix_canus,"
								"Pelomedusa_vereda),Cygnus_placidus)),Ciconia_mitratus),Nhandu_leucogaster),Spalax_thibetanus),"
								"(Teratolepis_resinifictrix,Trionyx_quadriocellata)),Limnodromus_completus),(Ameiva_marinus,"
								"Salvelinus_ochropus)),(Ctenotus_rosea,Iguana_bewickii)),(((Aegypius_vanellus,((Latastia_avocetta,"
								"Tropidurus_interiorata),Mesoplodon_leucopsis)),((((((Androctonus_caerulea,Felis_sepsoides),"
								"Streptopelia_alpinus),Sorex_cliffordii),Ovis_nasicus),Gyps_taxus),Spermophilus_oedicnemus)),"
								"(((((((((((Anolis_vegans,Pusa_caninus),Remiz_hungaricus),(Coturnix_agama,Psammophis_giganteus)),"
								"Vormela_monoceros),((Furcifer_subcinctus,Pterocles_squaterola),Macrorhamphus_mnemosyne)),"
								"(Canis_cocincinus,(Paramesotriton_serpentina,Vulpes_tuberculosus))),Rhinolophus_molurus),"
								"((((Enhydra_gibbosus,Morelia_nigriceps),Pseudemys_situla),Limnaeus_graeca),"
								"Vormela_exanthematicus)),((((((Callipogon_melonotis,(Lyrurus_martensi,Rhacophorus_terrestris)),"
								"((Chrysemys_guangxiensis,(Pelecanus_albirostris,Pterocles_atriceps)),(Fulica_geyri,Pelecanus_cyanea))),"
								"Pterocles_chinensis),((Gambelia_cygnus,Iguana_totanus),Gongylophis_ulikovskii)),Notophthalmus_melonotis),"
								"Felis_fimbriatus)),Gekko_taezanowskyi),Eurynorhynchus_auratus))),Phrynohyas_fuscatus),"
								"((((Almo_rosmarus,Eubalaena_casualis),Sorex_fusca),Apalone_nasicus),(Desnana_atra,Ephibolus_Bernicla)))"
								",Cyclemys_laevis),Falcipennis_fischeri),((Elseya_tarda,Halichoerus_cygnus),(Marmota_smaragdina,"
								"Ortigometra_erythropus))),Pratincola_coelestinus),((((Alauda_barroni,Hyperoodon_sp),"
								"(Norops_celer,Terpsihone_pholeter)),Paradoxornis_chukar),Haliaeetus_aegagrus)),"
								"((Emberiza_himantopus,(Oedura_crispus,Podoces_albatrus)),Zosterops_gecko)),Ctenosaura_caninus),"
								"((Emberiza_infrafrenata,Phrynosoma_flammea),Varanus_gecko))),(Bombyx_maihensis,Rhombomys_nigrolineatus))"
								",(((Buthacus_piscator,((Megaloperdix_bengalensis,Ptychozoon_ferox),Pyrrhocorax_sirtalis)),Panthera_mexicanum),"
								"Rissa_cinaedus)),(Ameiva_caudolineatus,(Cinclus_sphenocercus,(((Desnana_nupta,Eremophila_manul),Liasis_boa),"
								"Nyctaalus_standingii)))),((Dipsosaurus_hendricksoni,Philacte_bedriagai),Gallinago_caerulea)),Phrynops_prominanus),"
								"((((Bronchocela_bellii,Neolycaena_middendorffi),(((((Corytophanes_erythronota,Leiopython_comicus),"
								"(((Cuculus_flavescens,Pyrgilauda_torquatus),Xenochrophis_celer),(Cyclemys_tristis,Osteopilus_querquedula)))"
								",Haliaeetus_sibirica),(((Psammophis_sp,Theloderma_duplex),Varanus_hypoleucos),Ziphius_mongolica)),Nyroca_buccata))"
								",Seokia_tinctorius),Xenophrys_parahybana)),Chelus_placidus),(((((((((((((Budytes_docilis,"
								"(Poephagus_leschenaultii,Porzana_dauricus)),((((Diomedea_perdix,Glareola_livia),"
								"((Lampropeltis_novaeangliae,Notophthalmus_livia),Nyctaalus_paradisi)),Gypaetus_mystaceus),Testudo_cristata))"
								",Nemorhaedus_glacialis),Squaterola_pictus),Candoia_arcticus),Prunella_serricollis),Sternotherus_macularius)"
								",Pogona_iankowskii),((Euspiza_longipes,(Perdix_ocellatus,Rosalia_spaldingi)),Regulus_mexicana)),"
								"Latastia_monoceros),(Trachemys_deremensis,sibiricus_holbrooki)),Psammophis_cyanea),(Cygnus_wogura,Equus_ibis)))"
								",(((((((((((((((((Acanthosaura_cygnoides,(Basiliscus_deserti,Ctenotus_alpestris)),Lampropeltis_flavolineata)"
								",Rufibrenta_africanus),(Agama_tarandus,Terpsihone_zenobia)),((((((Alloporus_conicus,Phormictopus_erythronotus)"
								",((((Bombus_indicus,Pyrrhocorax_laticauda),Emydura_graeca),Streptopelia_dolosus),Terpsihone_docilis))"
								",Mogera_paradoxus),Limosa_belliana),Rhacophorus_taxus),Gonyosoma_margaritifera)),(Mabuya_bengalensis,"
								"Phoca_glacialis)),Dasypeltis_completus),Plegadis_glacialis),Aplopeltura_korschun),Rissa_aureola),"
								"(Chlidonias_nigrolineatus,Rhombomys_americanus)),Euspiza_haliaetus),((((((((Almo_communis,"
								"(Riparia_guentheri,Rosalia_belliana)),Atrophaneura_cygnoides),Parnassius_middendorffi),Sphenops_marmoratus)"
								",Vipera_chuatsi),Chalcides_bellii),(Buthus_cavirostris,Coleonyx_stimsoni)),(Homopholis_stylifer,"
								"Remiz_tridactylum))),(((((Alauda_albatrus,Osteopilus_leucocephala),(Lepidobatrachus_hilarii,Tadarida_papuana))"
								",((Heterodon_cynodon,Poephagus_albocinctus),((((Ketupa_czerskii,(Salamandra_leucocephala,Scorpio_mnemosyne))"
								",Ketupa_middendorffi),Rhacodactylus_serricollis),Sterna_mlokosiewiczi))),(Litoria_querquedula,Salmo_mehelyi)),"
								"((((((Archispirostreptus_brandtii,Tamias_cranwelli),Polypedates_wislizeni),(Chrysopelea_holbrooki,"
								"Circaetus_dives)),Bubulcus_bicoloratum),Micropalama_terrestris),Odonthurus_pallasii))),(Pandion_fluviatilis"
								",Phrynohyas_reinwardti)),sibiricus_porphyrio),Paramesotriton_mitratus)),Ninox_lesueurii),"
								"(((((((Bradypodion_brevipes,Tetrao_naumanni),Pterocles_versicolor),Budytes_catenifer),(Buteo_stagnalis"
								",Lamprolepis_cancerides)),((Branta_porphyrio,Dendrobates_keyzerlingii),Cardiocranius_jaculus)),"
								"((((Cottus_fischeri,Rhombomys_gibbosus),Ingerophrynus_colubrinus),Furcifer_ameiva),(Coturnix_scripta,"
								"(Megaptera_dactylisonans,Tetraogallus_cocincinus)))),Porphyrio_apollo)),Polypedates_hipposideros),"
								"Asthenodipsas_hendricksoni),Grammostola_getula),Net_reinwardti),Balaenoptera_rosea),(((((((Abantias_major"
								",Asthenodipsas_weberi),Dipsosaurus_longicollis),(((((((((((((((Abantias_musicus,(Alpes_leucoptera,"
								"((Aythya_leucoptera,Chlamydosaurus_crassicauda),Phormictopus_licin))),((Dafila_korschun,Pelomedusa_lobatus)"
								",Phelsuma_parvus)),Leuciscus_blythi),Leuciscus_guttifer),Scolopax_shadini),(((Eschrichtius_mitratus"
								",Leiocephalus_fragrans),Limnodromus_sudanensis),(Mochlus_vitulina,Psalmopoeus_leucophtalmos))),"
								"(Brachypelma_vastus,Mustela_galactonotus)),(Eurynorhynchus_caryocatactes,Leiopython_brongersmai))"
								",Monodon_cinaedus),((((Acanthoceros_nivalis,Geochelone_medici),Rangifer_schokari),Nemorhaedus_rosea)"
								",((((((((((Alectoris_lopatini,Siniperca_caucasicus),((Boiga_subrubrum,((Ceratophrys_standingii,"
								"(Chondropython_melonotis,(Hysterocrates_mongolica,Tringa_cancerides))),Rhynchophis_unicus)),Equus_gallicus))"
								",(((Bos_hemilasius,(Gavia_Anas,Kassina_scabra)),Chrttusia_garmani),((Lasiodora_aspera,Paramesotriton_floridana),Net_middendorffi))),Net_drapiezii),Bradyporus_comicus),Pyxicephalus_cygnus),Spalax_ibis),Enhydra_vertebralis),Ptychozoon_galeatus),Motacilla_interpres))),(Chrysopelea_monedula,Lepus_plumipes)),((((((Almo_jubata,Balaena_leucoptera),Tropidurus_heudei),Bubulcus_hendricksoni),((((Athene_riparia,Plethodon_radiata),Chondropython_semipalmatus),(Brachypelma_tricolor,(Sterna_rusticolus,Vipera_gemmicincta))),Glareola_aegyptia)),Crocodylus_leucostomum),(Latastia_angustirostris,Plegadis_coelestis))),((Alpes_arizonensis,Leiopython_schreibersi),(((Chamaeleo_bairdi,(Dryobates_limosa,(Lanius_cyanogaster,Otis_intermedia))),Phelsuma_vastus),Kaloula_turneri))),Brachypelma_martensi),(Homopholis_cristatellus,Motacilla_mugodjaricus))),(Net_dorsalis,Phoca_comicus)),(Acanthis_radiata,Pseudorca_placidus)),((Alectoris_brevipes,Chelydra_mackloti),(Buthacus_vittatus,Hysterocrates_melanoleucus))),(Homalopsis_hilarii,Rhombomys_dominicus))),(Citellus_aeruginosus,Sceloporus_geyri)),Iguana_coturnix),(((((((((Crotaphytus_rubida,Middendorffinaia_leucomystax),Tupinambus_tarda),(Fregilegus_hispida,Leiocephalus_longipennis)),(Himantopus_tarandus,Sus_enhydris)),Dyscophus_cranwelli),Pandion_brandtii),Sitta_monachus),Pseudemys_bicoloratum),Physignathus_septentrionalis)),Candoia_opimus),Zosterops_falcinellus),(Fuligula_minor,Rhombomys_rubida)),((((((((Balaena_ovata,(Chrysopelea_ruficollis,Mogera_enydris)),Tadarida_duplex),Phrynops_mutabilis),Melanocoryhpa_prasina),Fulica_ciliatus),Tadarida_resinifictrix),Cygnus_ochropus),Eschrichtius_hirundo)),(((((((((((((((((((Acanthosaura_colombianus,Gonyosoma_leschenaultii),Pyxicephalus_quadriocellata),((Pandinus_sanguinolentus,Sitta_leucomelas),Pterocles_guineti)),Nucifraga_galeatus),(Dryobates_scabra,((Emydura_crassidens,Paradoxornis_nigropalmatus),Uncia_nippon))),((((((((Antilope_indicus,(Dahurinaia_glareola,Rhinolophus_felderi)),(Desnana_scincus,Lasiodora_avosetta)),Ruticilla_sanguinolentus),(Phrynomerus_flammea,Python_tenuirostris)),Machetes_bicoloratum),Melanocoryhpa_garulla),Numenius_cambridgei),Tadorna_davidiana)),(Cervus_scrofa,((Hydrochelidon_zagrosensis,Melanocorypha_bengkuluensis),Ziphius_maldivarum))),((Ceratophrys_grunniens,(Litoria_galactonotus,(Pelomedusa_terrestris,Physignathus_brongersmai))),((Fuligula_caelebs,Upupa_fluviatilis),Lycaenopsis_vanellus))),Coregonus_hardwickii),(Cottus_ochropus,Eudramias_canus)),(Brachyramphus_garulla,(Enhudra_cioides,(Phyllopneuste_linaria,Rhinolophus_davidiana)))),Tadorna_clinatus),((((((((Ardea_tinctorius,Vanellus_himalayensis),Porzana_zenobia),Ethmostigmus_leiosoma),Osmoderma_flammea),Scincus_pygmaeus),(Fregilegus_asperum,Pelodiscus_ion)),Babycurus_fallax),Saga_caeruleus)),(Chlamydosaurus_kingii,(Chrysemys_godlewskii,Lepidobatrachus_rudicolis))),(((((Aegypius_multifasciata,((Bombina_triangulum,Natriciteres_stellatum),Chlamydotis_cygnoides)),Panthera_parvus),Pandion_medirostris),(Ardea_vegans,((Emberiza_anser,Otis_epops),Eurynorhynchus_lutra))),Net_acutus)),((Dipus_pyromelana,Rissa_pica),Enhydris_leiosoma)),(Coleonyx_musicus,Gallinago_barroni)),Dasypeltis_mexicanum),Turdus_nigropalmatus)),((((((Bombycilla_peregusna,((((Circus_vitticeps,((Lystrophis_altaica,Philomachus_vegans),(Procellaria_similis,Ursus_cinaedus))),Himantopus_cinerea),Eirenis_sujfunensis),Minipterus_dauricus)),Picus_medirostris),Haliaetus_caninus),((Oxyura_solitaria,Syrrhaptes_albopillosum),Phrynomerus_margaritifera)),(Hottentotta_pygargus,Riparia_colchicus)),Strepsilas_himalayensis)),(Acanthis_pulchra,(Coregonus_floridana,Eschrichtius_stylifer))),((Bufo_karelini,Rhesus_marmorata),Corallus_duplus)),((Eurynorhynchus_spinifera,Nipponia_mystaceus),Rufibrenta_insularis)),(((((Corvus_armeniacus,Motacilla_iguana),Xenophrys_azureus),Otis_grus),Ortigometra_helena),Pituophis_albigula)),Phormictopus_plumipes),Fuligula_comicus),Fulica_enhydris),Gonyosoma_orientalis),(Aquila_erythrogastra,Sphenurus_indica)),Atrophaneura_seemani),((Chrysopelea_brandtii,Cottus_turneri),Middendorffinaia_prominanus)),Salmo_mystaceus),((Babycurus_viscivorus,((Hyperoodon_similis,Podoces_bonasus),Turdus_citrsola)),((Castor_avocetta,Circus_eremita),Sceloporus_radiata))),Dahurinaia_coelestis),((((((((((((((((((((((((((((Acanthis_guineti,((Alopex_grunniens,(Chlamydotis_schreibersi,Rhacophorus_dendrophila)),Gonocephalus_wislizeni)),Alauda_chinensis),Bubulcus_bonasus),((((((Aegypius_triangulum,Megophrys_kuhli),Eublepharis_indica),(Gazella_bengalensis,Kinosternon_personata)),Liasis_hypoleucos),Rhodostethia_clypeata),((((Buthacus_cioides,Otocoris_caesius),Felis_placidus),(Charadrius_emarginatus,(Saxicola_leucomelas,Upupa_tadorna))),(Epipedobates_temminckii,Opheodrys_aristotelis)))),Gazella_melleri),(((Androctonus_maldivarum,(Petrocincla_capra,Tringa_pardalis)),Nyctaalus_indica),Iguana_agama)),Moschus_arvensis),((((((((Ahaetulla_diadema,Casarca_sanguinolentus),Norops_galeatus),(((Athene_minutus,Holodactylus_rhymnus),Panthera_semipalmatus),(Castor_solitaria,Parus_ammon))),Hirundo_heliaca),Saga_auriculatus),Gypaetus_collaris),(Ethmostigmus_mongolica,Phasianus_oedicnemus)),((Cyclagras_savignii,Lasiodora_longipes),(Mesoplodon_flammea,Plegadis_mexicanum)))),(Boa_molurus,Picus_tarda)),Archispirostreptus_terrestris),Sericinus_sphenocercus),((((((((((((((((((((((((((((((((((((((((((((Acanthis_mystacinus,((((((Actitis_dentata,Ceratophrys_subcinctus),((((((((((((((((Actitis_mitratus,Ophisops_noctua),((((((((Ambystoma_taezanowskyi,Dafila_resinifictrix),Pachydactylus_oedicnemus),Notophthalmus_tinnunculus),Hysterocrates_leucostomum),Hydrosaurus_perdix),(Neolycaena_jaspidea,Pituophis_Bernicla)),(Enhydra_pyromelana,Heteroscodra_gregaria)),(Gonyosoma_longipes,Uroplatus_insularis))),(((Calidris_caeruleus,(Circus_melleri,Neolycaena_kingii)),((Gerrhosaurus_acutus,Rhesus_maldivarum),Oedura_melanostictus)),(Capra_schreibersi,(Carabus_medirostris,(Eubalaena_caniceps,Falco_clypeata))))),Gazella_diffidens),(Apodora_corone,((Machetes_castaneus,Megophrys_caryocatactes),Neolycaena_caesius))),((Castor_temminskii,(((Elaphe_musculus,Paradoxornis_brandtii),Falcipennis_blythi),(Lystrophis_gregarius,Rhombomys_albigula))),Chlamydotis_situla)),(Bradyporus_ruficollis,((Chelus_tinnunculus,Sturnus_climacophora),Lycodon_grus))),Lepidobatrachus_arseniavi),(Chalcides_tridactylum,Ortigometra_capra)),((Babycurus_sp,(Heteroscodra_colombianus,(Paramesotriton_nigrolineatus,Telescopus_rubida))),Ovis_cherrug)),Pogona_proteus),(Lycaenopsis_nelsonii,(Spalerosophis_margaritifera,Terpsihone_auratus))),(((Chettussia_imperator,Leiocephalus_spaldingi),Phrynohyas_caesius),Coturnix_cachinans)),Megaloperdix_dominus),((Amphiuma_tolai,(Aphonopelma_graculus,Litoria_doriae)),(Epicrates_taezanowskyi,Fregilegus_undulata))),(((Alectoris_geniculata,Dyscophus_atra),Epipedobates_macularius),Hydrochelidon_guineti))),Phrynocephalus_dives),Tupinambus_nasuta),Apalone_prominanus),Xenophrys_sphenocercus)),Middendorffinaia_strepera),((((((Anser_bicoloratum,(Diomedea_fuscus,Minipterus_daurica)),((((((Apodora_aureola,Gerrhosaurus_glacialis),Recurvirostra_viridescens),((((Leiocephalus_cygnoides,Lycodon_himalayensis),(((Mustela_fulvus,Terpsihone_ferox),Phoca_melanoleucus),Ovis_martensi)),Siniperca_ignicapillus),Pareas_kuhli)),Cyriopagopus_sauritus),Chrysemys_bicinctores),Cervus_interiorata)),Phelsuma_durus),(((Lobipes_lopatini,Phrynomerus_glaucescens),Nucifraga_verrucosus),Vanellus_bicoloratum)),Thamnophis_rosea),Hyla_hassanica)),Synthliboramphus_filipjevi),(((Aegialifes_fuellebornii,Ahaetulla_carnivorus),((Basiliscus_veredus,Sphenops_himantopus),Enhydris_hypomelus)),(Citharacanthus_guttifer,(Rhacodactylus_imperator,Thecla_versicolor)))),Rhabdophis_cavimanus),((((((((((Alcedo_jaculus,Paraphysa_nyroca),Ingerophrynus_intermedia),(((Athene_bicoloratum,Trapelus_belliana),Eucratoscelus_scabra),Xenophrys_gigas)),Furcifer_gobio),Mogera_fluviatilis),(((Alloporus_fusca,(Chlamydosaurus_clinatus,Opheodrys_livia)),Recurvirostra_daurica),Eulabeia_tridactylum)),Charadrius_mnemosyne),Chlamydosaurus_bicinctores),Canis_amethistina),Scolopax_erythronotus)),Sternotherus_comicus),Uroplatus_undulata),(((((((((((((((((((Aegypius_smithii,Margaritifera_taezanowskyi),((Calidris_pygmeus,(Candoia_leschenaultii,Tropidurus_thibetanus)),((Callipogon_viridis,Ephibolus_arseniavi),Phyllopneuste_dennysii))),((Bombina_ovata,Cynops_livia),Micropalama_chuatsi)),Xenophrys_ussuriensis),((Anthropoidae_peregusna,Holodactylus_marinus),Trachemys_kraepelini)),Otis_martensi),Eubalaena_dorsalis),((((((((((((((((((((((Aix_stagnalis,(((Eirenis_leucomystax,Heteroscodra_mehelyi),(Leiocephalus_opimus,Rhesus_nyroca)),Lycodon_leucogaster)),Melanocorypha_guineti),Norops_crassidens),Citellus_major),Holaspis_ibera),Capra_gecko),(Archispirostreptus_saiga,Pterocles_conicus)),(((((Anthropoidae_quinquestriatus,Phormictopus_tristis),Sorex_quinquetaeniata),(((Chondropython_tentaculatum,Underwoodisaurus_garmani),((Enhydris_auratus,Rhynchaspis_gemmicincta),Tropidurus_miliaris)),Dendrobates_quinquetaeniata)),Clemmys_mitratus),(Hottentotta_ruficollis,Ketupa_ibera))),(Pterocles_alpinus,Pusa_plumipes)),((((((Anthropoides_iguana,Nyctixalus_jaculus),Tadorna_scabra),Canis_jaspidea),((Camptoloma_communis,Hirundo_mackloti),Hemiscorpius_quadrivirgata)),Net_angulifer),Leiocephalus_fimbriatus)),(Eutamias_erythronota,Gekko_chuatsi)),Kinixys_mlokosiewiczi),((((Anser_lineatus,((((Bombina_caniceps,Ceratophrys_variegatus),Monticola_riparia),(Myotis_alpinus,Uncia_caerulea)),(Lycodon_lopatini,Odonthurus_adspersus))),Furcifer_euptilura),Pyrgilauda_hypoleucos),Lamprophis_sauromates)),((((((Aplopeltura_guttifer,((((Chlidonias_aegagrus,Dendrobates_platyrhinos),(Chlidonias_heudei,Erpeton_caelebs)),Vormela_salvator),Kinosternon_subminiatus)),Coleonyx_anatina),Epicrates_situla),Basiliscus_zagrosensis),((((Chelodina_squamatus,Strepsilas_insignis),(Laudakia_cavimanus,Ziphius_erythronota)),((((Hydrochelidon_fusca,Leptobrachium_mykiss),Recurvirostra_ion),((Milvus_verrucosus,Uroplatus_prominanus),Polypedates_cioides)),Sternotherus_obsoleta)),Lobipes_durus)),Regulus_enydris)),Callipogon_ovata),Oligodon_reinwardti),(((((((((((((Alauda_zagrosensis,((Aythia_squamatus,Chrysemys_irregularis),Lyrurus_rubida)),((((((Avicularia_agama,(Damon_macqueni,Lanius_jubata)),(Lanius_tricolor,Syrrhaptes_percnopterus)),Pelodytes_brachydactyla),((Cyclagras_cinclus,(((Cygnus_dispar,((Emberiza_melleri,Sphenops_brongersmai),Xenopeltis_mehelyi)),Gypaetus_exquisita),Gekko_sanguinolentus)),Saiga_fluviatilis)),Falco_ruficollis),Bubulcus_ocellatus)),(((((((((((Allactaga_leucogaster,((Anolis_exanthematicus,Parnassius_ovata),Spalax_docilis)),Scincus_maurus),Enhydra_caryocatactes),Ketupa_wogura),Avicularia_atrigularis),Eubalaena_sudanensis),Underwoodisaurus_moschata),(Gypaetus_battersbyi,(Rangifer_bonasus,Xenophrys_grunniens))),((Avicularia_argentatus,Phrynops_gallicus),(Buthacus_caelebs,Eremophila_dubius))),(Columba_guttifer,Ortigometra_quinquestriatus)),(Circaetus_personata,(Rhacophorus_longicollis,Rhesus_quinquestriatus)))),Oligodon_albigula),Cygnopsis_vitulina),((Calidris_hilarii,(((Dahurinaia_himantopus,Limosa_albicilla),(Rosalia_gallinago,Uromastyx_ferruginea)),Pagophila_armeniacus)),Pelodiscus_cioides)),(((((Anas_tetrax,Philomachus_calvus),((((Bubulcus_altaicus,Crotaphytus_pygargus),(Sturnus_gigas,Ziphius_communis)),Homalopsis_Jankowskii),Tringa_mirabilis)),Middendorffinaia_cinaedus),(Mustela_taxispilota,Psalmopoeus_carbo)),(((Bubulcus_walti,(Madagascarophis_elaphus,Ruticilla_belliana)),Tropidurus_hosii),Lamprophis_means))),(Bronchocela_dorsalis,((Cynops_insignis,(Gazella_isabellina,Strepsilas_fissipes)),Saxicola_avicularia))),((Callipogon_scincoides,((Chlamydosaurus_fallax,Rhinolophus_tigrinus),Phrynomerus_ovata)),Fulica_canagica)),Rufibrenta_battersbyi),Lutra_ocellatus),Margaritifera_fulvus),Enhydra_caudata)),Bombyx_ignicapillus),Candoia_pyromelana),(Basiliscus_mitratus,(Latastia_elaphus,Moschus_pulcher))),(Gongylophis_japonensis,Ursus_cyanogaster)),Spermophilus_holbrooki)),Siniperca_alba),Turdus_eburnea),Rhamphiophis_agama),(((((((((Ahaetulla_flavigularis,(Phalacrocorax_ammon,Pterinochilus_stejnegeri)),Eubalaena_Anas),((Homalopsis_leucophyllata,Phrynohyas_venulosa),(Mabuya_rubida,Nerodia_cepediana))),Cardiocranius_atthis),Regulus_gregarius),Hemitheconyx_plathyrhychos),Euspiza_indicus),Uncia_monorhis),(Avicularia_mackloti,Moschus_flavolineata))),Hydrochelidon_gemmicincta),Ciconia_leucotus),Sorex_bairdii),(Alopex_proteus,(((Callipogon_hemilasius,Macrorhamphus_americanus),Dipsosaurus_blakistoni),Chlamydotis_pulchripes))),(Chelydra_ciliatus,Spizaetus_blythi)),Paraphysa_veredus),(((Boa_ciliatus,Ketupa_turneri),((Bufo_piscator,(((Charadrius_bengalensis,Epicrates_vegans),Lasiodora_teguixin),Gavia_macrops)),Ortigometra_iankowskii)),Scincus_albatrus))),(Falcipennis_armeniacus,Holodactylus_fallax)),Tetrao_baeri),Tamias_maritimus),((Kaloula_aureostriata,Spizaetus_aceras),Thecla_atthis)),Pyrrhocorax_sieboldii),((Accipiter_teniotis,Branta_lopatini),((((((((((((((((((Aegypius_gebleri,((((((((((((Chrysemys_major,(Geochelone_nigra,Thecla_castus)),(Gavia_lutris,Vulpanser_australis)),Crotaphytus_americanus),((Fulica_striatus,(Gekko_elaphus,Lasiodora_belliana)),Ketupa_unicolor)),Petrocincla_exquisita),Elseya_tatarica),Heterodon_pygmeus),Nerodia_rubicola),Phylloscopus_arseniavi),(((((Coenobita_erythrogastra,Gazella_scrofa),Hottentotta_lobatus),Erpeton_cingulata),Ninox_pyromelana),Philothamnus_cristatellus)),Nemorhaedus_dendrophila),Minipterus_bengkuluensis)),Moschus_helena),Gonyosoma_clypeata),Buteo_leucotus),(((((((Alaus_tricolor,Litoria_auriculatus),(((((((Anser_nebularia,Zosterops_placidus),Ctenosaura_wumuzusume),Morelia_subglobosa),Ziphius_fuscatus),Antilope_reticulatus),(Bombus_paradisi,Squaterola_subglobosa)),Bombina_mutabilis)),(((Alectoris_leucopsis,Elaphe_mykiss),Perdix_clypeata),(((Canis_proteus,Sceloporus_gordoni),Nhandu_clarus),(Hadogenes_wislizeni,(Psammophis_turtor,Rhynchophis_auratus))))),Epipedobates_ridibundus),(Parnassius_orientalis,Recurvirostra_sanguinolentus)),Bombycilla_subglobosa),(Netta_boyciana,Tamias_nigriceps))),(Homalopsis_graeca,((Nyctixalus_spinifera,Regulus_baibacina),Saga_chrysargos))),(Canis_lavaretus,Onychodactylus_carnifex)),(((((Capreolus_collaris,((Ctenotus_capra,(Euspiza_serricollis,Rhacodactylus_marmorata)),Mabuya_mexicanum)),Procellaria_caudata),Teratolepis_flavescens),Eublepharis_caudata),Hadogenes_hirundo)),(Candoia_nivicola,Citellus_kraepelini)),Testudo_aleutica),(Oceanodroma_chamaeleontinus,Tamias_molurus)),((((((Alopex_crispus,Neolycaena_zenobia),Nyctixalus_rufinus),((Eudrornias_semipalmatus,Siniperca_smithii),Scaphiopus_alpina)),Hydrosaurus_flavomaculatus),(((Archispirostreptus_vitulina,Gypaetus_mitratus),Kaloula_marcianus),Milvus_cancerides)),Corytophanes_ridibundus)),Vanellus_drapiezii),Vipera_piscator),((Mylopharyngodon_leucomystax,Vanellus_terrestris),Parnassius_clinatus)),Nucifraga_bobac),(Citellus_atra,(Euspiza_plumifrons,Mareca_hungaricus))),Buthacus_moschiferus))),Lepus_pardalis),(Sorex_karelini,Sphenurus_cianeus)),Numenius_merganser),Acanthosaura_sauromates),((Actitis_godlewskii,Ortigometra_cepediana),(Lamprolepis_fragrans,Lepus_vittatus))),((((((((Ardea_weliczkowskii,Ctenosaura_anser),((((Epipedobates_rosea,Minipterus_armeniacus),Hyla_cepediana),(Iguana_avocetta,Leiopython_fimbriatus)),Phrynomerus_heliaca)),Cygnus_boulengeri),(Elseya_ferruginea,Gavia_torquatus)),Monodon_pica),Middendorffinaia_carnifex),Citellus_pholeter),Canis_fuliginosus)),Myotis_relictus),(((((((Arenaria_mackloti,Monticola_leucogeranus),Calotes_glottis),(Calidris_sibirica,Hemiscorpius_vitulina)),Sorex_castaneus),(Mergus_guttifer,Phasianus_rufus)),Macrorhamphus_hendricksoni),((Bradyporus_leucomystax,Cyclemys_sphenocercus),Lycodon_africanus))),((Agama_dominus,(Holaspis_semipalmatus,Vanellus_labiatus)),(((Alaus_trianguligerus,Enhudra_melonotis),Diomedea_caelebs),Tadarida_fuscus))),Net_maurus),Synthliboramphus_ibis),Rhamphiophis_melanostictus),((Kinixys_marmoratus,Ptyodactylus_carolinensis),Physignathus_veredus)),((((((((((Alectoris_laticauda,Rhesus_gregaria),Pachytriton_albertisii),Hottentotta_parreyssi),(Bufo_lopatini,Sitta_quinquestriatus)),Odonthurus_longipennis),Sericinus_falcipennis),Pleurodeles_diffidens),(((Bos_gallinago,(Eremophila_flavigularis,(Hydrochelidon_wislizeni,Oedura_relictus))),Bubulcus_nipalensis),Middendorffinaia_kingii)),Trionyx_sepsoides),Saga_cristatus)),Ethmostigmus_mehelyi),Ruticilla_nigropalmatus),(((((((((((((Acanthogonatus_savignii,(((Capreolus_monacha,(Mustela_tridactylum,Nyctixalus_interiorata)),(Eudramias_leucomelas,Phoca_alcinous)),Mochlus_rufodorsata)),Pachydactylus_molurus),(Phoca_hemilasius,Pyxicephalus_aristotelis)),Osmoderma_barroni),Brachypelma_zagrosensis),Cygnus_blythi),Lasiodora_classicus),Saga_pedo),Psammophis_longipennis),Hyla_lesueurii),Bombyx_bengkuluensis),Tringa_sujfunensis),(Alloporus_pelagicus,Hemitheconyx_trigonopodus))),(((((((((((((((Acanthoceros_caninus,(((((Alcedo_hypomelus,Felis_avinivi),Spalerosophis_aegyptia),Vanellus_spaldingi),Heteroscodra_vipio),((((((Carabus_subglobosa,Lanius_vanellus),Mogera_bobac),Saiga_amboinensis),Corytophanes_acutus),(Ingerophrynus_hendricksoni,(Squaterola_gratiosa,Tiliqua_pugatshuki))),Oxyura_turneri))),Acheron_gordoni),((Aegialifes_murinus,Alpes_hyemalis),Odonthurus_grupus)),((Acipenser_javanica,Cygnopsis_tetrax),(Dendrelaphis_brandtii,Oenanthe_maculata))),Cyclagras_schrencki),(Cuculus_solitaria,(Pelusios_multituberculatus,Platalea_tricolor))),Philothamnus_cavirostris),Athene_subcinctus),Hydrochelidon_cherrug),Eulabeia_similis),Megaptera_truncatus),((((((Alauda_fusca,Gongylophis_ibis),Apalone_albirostris),Pica_sauromates),Python_jacksoni),(Gallinago_smithi,Holaspis_gregaria)),(Kassina_molurus,Philomachus_physalus))),(((((((((((((((((((((((((((((((((Acanthogonatus_leucogaster,Chrttusia_cranwelli),((((((((((Ameiva_irregularis,Lamprophis_griseus),Odobenus_niloticus),(Argynnis_erythropus,Dendrobates_nivicola)),(Callipogon_flava,(Coturnix_ibis,Phrynohyas_caelebs))),((Bradypodion_leucophyllata,(Pyxicephalus_atrigularis,Sphenurus_citrsola)),Triturus_serricollis)),Nucifraga_venulosa),Asthenodipsas_occitanus),Odonthurus_unicolor),Marmota_pulchra),Mergus_platycephala)),((((((((((((((Acanthogonatus_lopatini,(((Chamaeleo_leucophyllata,Perdix_glottis),(Dipus_maldivarum,Lutra_maculata)),Tadarida_clericalis)),Vipera_keyzerlingii),Micropalama_caelebs),Platalea_dahurica),Enhudra_euptilura),(((Acanthoscurria_exquisita,(Anas_aeruginosus,Scolopendra_homeana)),(Dipsosaurus_caelebs,Paramesotriton_dione)),(Ephibolus_castaneus,Ophisops_soloensis))),((((Ambystoma_medici,Corallus_docilis),Testudo_bedriagai),Chamaeleo_dolosus),Dipus_lineatus)),(Gongylophis_tentaculatum,Minipterus_ussuriensis)),Rhesus_grandis),(Capella_standingii,Gongylophis_calligaster)),Perdix_rupestris),Hyperoodon_kuhli),(((((((((((((Accipiter_krueperi,Poephagus_piceus),((((Caiman_azureus,Trapelus_epops),Spalerosophis_campestris),Cuon_deremensis),Emberiza_chuatsi)),Bradyporus_cancerides),((((((((((((((Carabus_monachus,((Lyrurus_communis,Ruticilla_septentrionalis),sibiricus_leptochelis)),Salamandra_paradoxus),(((((Cyriopagopus_leucocephala,Trionyx_alpestris),Milvus_notaeus),Mabuya_madagascariensis),Telescopus_taxispilota),Holaspis_franckii)),((((Chlamydotis_rosmarus,Leptopelis_olivacea),Dendrobates_vereda),(Hottentotta_carnivorus,Leptopelis_leporosum)),((Dryobates_yeltoniensis,(Holodactylus_parreyssi,Phrynops_maculata)),Spizaetus_maculata))),Heterodon_soloensis),Zosterops_shadini),Madagascarophis_elegans),Mesoplodon_minor),Leiurus_aegyptia),Hottentotta_gobio),Hyla_ulikovskii),(Chrttusia_isabellina,Damon_alpestris)),Querquedula_deremensis),Eucratoscelus_kuhli)),(Dafila_cioides,Desnana_schokari)),Columba_unicus),(Anthropoidae_cenchria,(((Colaeus_alterna,Mabuya_fragrans),Homopholis_purpurascens),(Passer_crucigera,Triturus_crassidens)))),Philomachus_yeltoniensis),(Bradypodion_fiber,Calotes_iguana)),Scorpio_tadorna),Cyclemys_erythronota),Ambystoma_bedriagai),(((((((Alectoris_gemmicincta,Colaeus_leucostomum),Lepus_mycterizans),Pareas_helena),Charadrius_femoralis),Ardea_eremita),Nyctixalus_hardwickii),Chamaeleo_nasuta))),(Aix_guttata,Philothamnus_hasselquistii))),Pachydactylus_guentheri),(Columba_govinda,(Elseya_erythrogastra,Lagenorhynchus_leucomystax))),Melanocorypha_amboinensis),Anodonta_mongolica),Gazella_rubida),Eulabeia_triangulum),(Acheron_calvus,(Odobenus_caerulea,Phrynocephalus_kuhli))),Spalax_platyrhinos),Parnassius_leucogaster),Cypselus_fissipes),((Anthropoides_penelope,Haliaeetus_aureostriata),Pandinus_Jankowskii)),Hyla_squamatus),Hadrurus_comicus),Tryngites_rutilans),Cuon_paradisi),(Emberiza_pugnax,Mogera_manul)),((Cypselus_garmani,Pyxicephalus_breitensteini),(Parnassius_filipjevi,Vipera_gigas))),Pelomedusa_guineti),Anser_nipalensis),((((Boiga_subruficollis,((Hydrochelidon_serricollis,Triturus_pica),(Opheodrys_ion,Pituophis_carbonaria))),Branta_ceterus),(Fuligula_americanus,Oxyura_laticauda)),Lamprolepis_squaterola)),Ctenosaura_albigula),Odonthurus_fasciolata),Agama_anser),(Bradyporus_plathyrhychos,(Cygnus_Anas,Sceloporus_hendricksoni))),(((((((((Ambystoma_novaeguineae,(((((((Argynnis_prominanus,Lamprolepis_xanthocheilus),(Budytes_glacialis,Charadrius_penelope)),(Certhia_subniger,Furcifer_citreola)),Megaptera_paradisi),(Lobipes_adamsii,Mogera_scincoides)),Vanellus_korschun),Mochlus_stimsoni)),((Gecarcinus_leucoryphus,((Parabuthus_purpurascens,Synthliboramphus_hongkongensis),Tupinambus_lavaretus)),Plegadis_lutra)),Plegadis_enydris),((Argynnis_calidris,(Lagenorhynchus_albicilla,(Saga_quadriocellata,Terpsihone_zagrosensis))),Melanocoryhpa_docilis)),(((((((((Bronchocela_mugodjaricus,Scorpio_hemilasius),Hydrosaurus_eximia),Spizaetus_longipes),Strepsilas_leucoryphus),(Eubalaena_truncatus,Mergus_rufina)),Chelydra_aristotelis),Vanellus_nivicola),Cuon_geyri),Nemorhaedus_limosa)),(Callipogon_bewickii,Totanus_schreibersi)),Cuora_caudicinctus),Mogera_clericalis),Larus_grandis)),Rosalia_dione),Candoia_armeniacus),Liasis_caelebs),Thymallus_gregarius),Anolis_smithii)),Anthropoidae_epops),Norops_hemilasius)),Teratoscincus_hyperboreus),Buthus_scalaris),((((((Chelus_nebrius,Lutra_tentaculatum),Lutra_tolai),Grampus_lobatus),Dipus_japonica),Columba_lutris),(Rhacophorus_uluguruensis,Teratolepis_dione))),(Apalone_quinquestriatus,Cypselus_brachydactyla)),(Eucratoscelus_nigropalmatus,(Seokia_means,Sorex_orientalis))),Hadrurus_gemmicincta),((((((((((Acanthosaura_ussuriensis,((Athene_turneri,((Egretta_edulis,Pelodiscus_caudatus),Lycodon_coturnix)),Machetes_salvator)),Aegypius_ornata),Egretta_salamandra),(((Aythia_stimsoni,(((Lagenorhynchus_carbonaria,Megaloperdix_ruficollis),Tetrao_peregusna),Parus_variegatus)),Motacilla_campestris),(Hydrochelidon_caninus,Salamandra_cingulata))),((((((((Aegypius_falcinellus,((((((((Buthacus_taezanowskyi,(((Chrysemys_walti,Pandion_helvetica),(Phrynosoma_tinctorius,Sphenops_javanica)),Totanus_graeca)),Liasis_flavolineata),Geochelone_ferrumequinum),Chrysopelea_thibetanus),Osmoderma_constricticollis),Castor_pusilla),Rissa_flavescens),Philacte_ferox)),Epicrates_flava),Dendrelaphis_leucocephala),Emydura_gemmicincta),Sternotherus_sebae),(((((Athene_anatina,Selenocosmia_castaneus),Mesoplodon_naumanni),Xenochrophis_maculatum),Enhydris_ladogensis),(Hydrochelidon_schneideri,Scaphiopus_peregusna))),((((((Aythia_dauricus,Tiliqua_diffidens),(((Bos_arenarius,((Ingerophrynus_grandis,Middendorffinaia_griseus),Turdus_stagnalis)),Motacilla_cavimanus),(Pandion_ameiva,Theloderma_leucogeranus))),Kaloula_alpinus),Megaptera_uluguruensis),Bombus_terrestris),(Casarca_gemmicincta,((Eirenis_sp,Lepidobatrachus_heudei),Physignathus_infrafrenata)))),Spalax_guttifer)),Leiopython_subniger),(Rhynchophis_pulcher,sibiricus_alpestris)),Apodora_kuhli),(Apus_mykiss,Eudrornias_sphenocercus)),Saga_parvus)),Ophisops_variegatus),Circaetus_grossmani),((Megaloperdix_bedriagai,Syrrhaptes_aureostriata),(Pandion_vittatus,Pleurodeles_totanus)))),Numenius_isabellina),((((((Allobates_major,Tringa_arseniavi),(Hysterocrates_pictus,Otocoris_dactylisonans)),(((((Burhinus_totanus,Psalmopoeus_taxus),(((Chelydra_nippon,Osmoderma_livia),Nyctixalus_veredus),Scaphiophryne_lutra)),Motacilla_major),Tadarida_glottis),Lamprolepis_sebae)),((((((Calotes_tolai,Porphyrio_bairdi),Heteroscodra_sagrei),Gecarcinus_chrysaetus),Remiz_calligaster),Pseudemys_cliffordii),Phrynohyas_campestris)),(Leuciscus_wumuzusume,Platemys_thibetanus)),(Aplopeltura_scincus,Procellaria_amethistina))),Podoces_elaphus),Alaus_flavescens),Leiurus_maldivarum),(Chamaeleo_percnopterus,Colaeus_anachoreta)),Columba_gigas),Hadrurus_squaterola),(((((((((((Bos_cherrug,Spermophilus_armeniacus),Spizaetus_musicus),Megaptera_buccata),Phasianus_viscivorus),Gypaetus_rostratus),((Bos_crocodilus,((Chen_filipjevi,(Chrysemys_brongersmai,Hadogenes_collaris)),((((Enhydris_semipalmatus,Psammophis_musculus),Gecarcinus_stylifer),Sphenops_gallicus),Fulica_temminckii))),Sterna_niloticus)),Vipera_ulikovskii),((Limosa_milii,Plegadis_paradisi),Psalmopoeus_cristatella)),Leuciscus_fasciolata),Grammostola_longipes),((Bronchocela_savignii,(((Cuora_rufus,(Hemitheconyx_leucotus,(Himantopus_wogura,Lamprolepis_sieboldii))),(Gekko_calamita,Motacilla_gordoni)),Neophron_calamita)),Rhamphiophis_vegans))),Rhynchaspis_calamita),Holodactylus_squamatus),Chlidonias_torquatus),((((((((((Allactaga_medici,Leiurus_eburnea),Geochelone_leucocephala),Dendrelaphis_viridescens),Sericinus_guttifer),((Himantopus_subglobosa,Kinixys_schreibersi),Recurvirostra_cingulata)),(((((((Anodonta_grunniens,Rissa_ion),(Athene_stellio,(Cypselus_cyanus,Eulabeia_calligaster))),Monodon_fiber),((Gekko_difficilis,Rhombomys_anser),(Scolopax_caniceps,Sturnus_bonasus))),Tadorna_mitratus),Haplopelma_vulgaris),(Canis_meles,Enhudra_albopillosum))),Ophisops_clinatus),Spermophilus_montela),(Dendrobates_nyroca,Dipus_hasselquistii)),(Moschus_plathyrhychos,sibiricus_salei))),((Bradypodion_leucotus,Gavia_platycephala),Lycodon_climacophora)),Fuligula_cranwelli),Neophron_emarginatus)),Abantias_leuconotus);");
	}
	
 }