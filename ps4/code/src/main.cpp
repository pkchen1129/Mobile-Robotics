#include <2D_Graph_SLAM.h>
void Batch2D(){
    //readvertex ( indexed pose , initial )
    std::pair<vector<boost::optional<IndexedPose>>, Values::shared_ptr> pair_vertex = readVertex();
    //readEdge(indexed pose, factor())
    std::pair < vector<boost::optional<IndexedEdge>> ,vector<NonlinearFactor::shared_ptr> > pair_edge = readEdge();
    vector<boost::optional<IndexedEdge>> vec2 = pair_edge.first;
    vector<NonlinearFactor::shared_ptr> vec_factor = pair_edge.second;
    // Ways to traverse through vec and vec2;
    for (auto it = pair_vertex.first.begin() ; it < pair_vertex.first.end() ; it++){
        auto index_pose = *it;
        Key id = index_pose->first;
        Pose2 pose = index_pose->second;
        // cout << pose.theta() << endl;
    }
    for (auto it2 = vec2.begin() ; it2 < vec2.end() ; it2++){
        auto between_pose = *it2;
        Key id1 = between_pose->first.first;
        Key id2 = between_pose->first.second;
        // cout << id1 << ", " << id2 <<endl;
        Pose2 p = between_pose->second;
        // cout <<  << endl;
    }
    // Creating graph for optimization
    NonlinearFactorGraph::shared_ptr graph(new NonlinearFactorGraph);
    
    for (auto it3 = vec_factor.begin() ; it3 < vec_factor.end() ; it3++){
        NonlinearFactor::shared_ptr factor = *it3;
        graph->push_back(factor);// Pushing factor to a graph
    }
    NonlinearFactorGraph initial_graph_WithPrior = *graph;

    // Adding initial Prior // 
    noiseModel::Diagonal::shared_ptr priorModel = noiseModel::Diagonal::Variances(Vector3(1e-6, 1e-6, 1e-8));
    initial_graph_WithPrior.add(PriorFactor<Pose2>(0, Pose2(), priorModel));
    std::cout << "Adding prior on pose 0 " << std::endl;

    // Optimization //
    int maxIterations = 100; // default
    GaussNewtonParams params;
    params.setVerbosity("TERMINATION"); // parameters
    std::cout << "Optimizing the factor graph" << std::endl;
    params.maxIterations = maxIterations;
    Values::shared_ptr initial = pair_vertex.second; // Putting initial into 
    GaussNewtonOptimizer optimizer(initial_graph_WithPrior, *initial, params);
    Values result = optimizer.optimize();
    std::cout << "Optimization complete" << std::endl;

    // output file
    string outfile = "../dataset/out_batch_2D.txt";
    writeg2o(initial_graph_WithPrior, result, outfile);
}
void Incremental2D(){
    //readvertex ( indexed pose , initial )
    std::pair<vector<boost::optional<IndexedPose>>, Values::shared_ptr> pair_vertex = readVertex();
    Values::shared_ptr initial = pair_vertex.second; // Putting initial into 
    //readEdge(indexed pose, factor
    std::pair < vector<boost::optional<IndexedEdge>> ,vector<NonlinearFactor::shared_ptr> > pair_edge = readEdge();
    vector<boost::optional<IndexedEdge>> vec2 = pair_edge.first;        
    vector<NonlinearFactor::shared_ptr> vec_factor = pair_edge.second;
    
    
    // Create an iSAM2 object. Perform partial relinearization/recording at each step.
    ISAM2Params parameters;
    parameters.relinearizeThreshold = 0.01;
    parameters.relinearizeSkip = 1;
    ISAM2 isam(parameters);

    // Create a Factor Graph and Values to hold the new data
    NonlinearFactorGraph graph;
    // NonlinearFactorGraph::shared_ptr graph(new NonlinearFactorGraph);
    
    // Add the prior
    noiseModel::Diagonal::shared_ptr priorModel = noiseModel::Diagonal::Variances(Vector3(1e-6, 1e-6, 1e-8));
    graph.add(PriorFactor<Pose2>(0, Pose2(), priorModel));
    std::cout << "Adding prior on pose 0 " << std::endl;
    // Add the init
    Values initialEstimate; // Create values variable to place poses
    initialEstimate.insert(pair_vertex.first[0]->first, pair_vertex.first[0]->second); //Insert pose 0 from g2o file
    
    Values result;
    // update ISAM2 with the NonlinearFactorGraph and Values variables
    isam.update(graph, initialEstimate);
    result = isam.calculateEstimate(); // What format would this be?

    for(int i = 1 ; i < pair_vertex.first.size() ; i++){
        // Clear graph and init value
        graph.resize(0);
        initialEstimate.clear();

        // Add the current node to init value using result.at(i) + odometry(should transform to global matrix)
        initialEstimate.insert(pair_vertex.first[i]->first, pair_vertex.first[i]->second);

        // noiseModel::Diagonal::shared_ptr model = noiseModel::Diagonal::Sigmas(Vector3(0.2, 0.2, 0.1));
        // auto l1Xl2 = gtsam::DerivedValue::between(pair_vertex.first[i]->second, pair_vertex.first[i-1]->second);
        // graph.emplace_shared<BetweenFactor<Pose2> >(i-1, i, Pose2(2, 0, M_PI_2), model);
        // initialEstimate.insert(i+1, initialEstimate.at<Pose2>(i) * l1Xl2);

        for(int j = 0 ; j < pair_edge.first.size() ; j++){ // 
            // find edge that connects to the current nodes and previous added nodes
            // cout << pair_edge.first[j]->first.second << ", ";
                // cout << "hi" << endl;
            if(pair_edge.first[j]->first.second == i){ // id2 == current node
                // add the edge to the graph
                graph.push_back(pair_edge.second[j]); // Graph push_back factors
            }
        }

        // isam update
        isam.update(graph, initialEstimate); // It can't update.
        // isam calculate estimate
        result = isam.calculateEstimate(); // What format would this be?
        
    }
    // output file
    std::cout << "Output File... " << std::endl;
    string outfile = "../dataset/out_isam2_2D.txt";
    writeg2o(graph, result, outfile);
    
}
void Batch3D(){
    string filename = "../dataset/parking-garage.txt";
    const auto factors = p3DFactors(filename);
    NonlinearFactorGraph::shared_ptr graph(new NonlinearFactorGraph);
    for (const auto& factor : factors) {
        graph->push_back(factor);
    }

    const auto poses = p3DPoses(filename);
    Values::shared_ptr initial(new Values);
    for (const auto& key_pose : poses) {
        initial->insert(key_pose.first, key_pose.second);
    }

    NonlinearFactorGraph graphWithPrior = *graph;
    // Adding initial Prior // 
    noiseModel::Diagonal::shared_ptr priorModel = noiseModel::Diagonal::Variances((Vector(6) << 1e-6, 1e-6, 1e-6, 1e-4, 1e-4, 1e-4).finished());
    Key firstKey = 0;
    for(const Values::ConstKeyValuePair& key_value: *initial) {
        std::cout << "Adding prior to g2o file " << std::endl;
        firstKey = key_value.key;
        graphWithPrior.add(PriorFactor<Pose3>(firstKey, Pose3(), priorModel));
        break;
    }


    // Optimization //
    int maxIterations = 100; // default
    GaussNewtonParams params;
    params.setVerbosity("TERMINATION"); // parameters
    std::cout << "Optimizing the factor graph" << std::endl;
    GaussNewtonOptimizer optimizer(graphWithPrior, *initial, params);
    Values result = optimizer.optimize();
    std::cout << "Optimization complete" << std::endl;

    // output file
    string outfile = "../dataset/out_batch_3D.txt";
    writeg2o(graphWithPrior, result, outfile);
}
void Incremental3D(){
    // 3D Incremental    
    string filename = "../dataset/parking-garage.txt";
    const auto factors = p3DFactors(filename);
    // Iterate method:
    // factors[i]->key1()
    // factors[i]->key2()

    const auto poses = p3DPoses(filename);
    // Iterate Method: 
    // poses.find(i)->first
    // poses.find(i)->second // Quaternion: x,y,z + (qw, qx, qy, qz)


    // Create an iSAM2 object. Perform partial relinearization/recording at each step.
    ISAM2Params parameters;
    parameters.relinearizeThreshold = 0.01;
    parameters.relinearizeSkip = 1;
    ISAM2 isam(parameters);

    // Create a Factor Graph and Values to hold the new data
    NonlinearFactorGraph graph;
    
    // Add the prior
    noiseModel::Diagonal::shared_ptr priorModel = noiseModel::Diagonal::Variances((Vector(6) << 1e-6, 1e-6, 1e-6, 1e-4, 1e-4, 1e-4).finished());
    graph.add(PriorFactor<Pose3>(0, Pose3(), priorModel));
    std::cout << "Adding prior on pose 0 " << std::endl;
    
    // Add the init
    Values initialEstimate; // Create values variable to place poses
    initialEstimate.insert(poses.find(0)->first, poses.find(0)->second); //Insert pose 0 from g2o file
    Values result;
    
    // update ISAM2 with the NonlinearFactorGraph and Values variables
    isam.update(graph, initialEstimate);
    result = isam.calculateEstimate(); // What format would this be?

    for(int i = 1 ; i < poses.size() ; i++){
        // Clear graph and init value
        graph.resize(0);
        initialEstimate.clear();

        // Add the current node to init value using result.at(i) + odometry(should transform to global matrix)
        initialEstimate.insert(poses.find(i)->first, poses.find(i)->second);
        noiseModel::Diagonal::shared_ptr priorModel = noiseModel::Diagonal::Variances((Vector(6) << 1e-6, 1e-6, 1e-6, 1e-4, 1e-4, 1e-4).finished());
        for(int j = 0 ; j < factors.size() ; j++){ // 
            // find edge that connects to the current nodes and previous added nodes
            // cout << pair_edge.first[j]->first.second << ", ";
            if(factors[j]->key2() == i){ // id2 == current node
                // add the edge to the graph
                graph.add(factors[j]);
            }
        }

        // isam update
        isam.update(graph, initialEstimate); // It can't update.
        // isam calculate estimate
        result = isam.calculateEstimate(); // What format would this be?
    }
    // output file
    std::cout << "Output File... " << std::endl;
    string outfile = "../dataset/out_isam2_3D.txt";
    writeg2o(graph, result, outfile);
}
int main() {
    Batch2D();
    // Incremental2D();
    Batch3D();
    Incremental3D();



    return 0;

}


