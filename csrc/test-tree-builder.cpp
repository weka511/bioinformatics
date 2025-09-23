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
	
	 SECTION("Test parser with labels") {
			shared_ptr<Parser::Tree> tree = parser.parse("(A,B,(C,D));");
			tree->descend(builder);
			shared_ptr<Node> new_tree = builder->get_result();
			cout << *new_tree << endl;
			new_tree->visit(displayer);
			// cout<<displayer->get_string()<<endl;
	 }
	 
	 SECTION("Test parser with labels and lengths(1)") {
		shared_ptr<Parser::Tree> tree = parser.parse("(A:2.1,B:1.2,(C:3.0,D:0.2));");
		tree->descend(builder);
		// cout<<builder->get_string()<<endl;
	}
}