

/**
 * @file 2D_Graph_SLAM.cpp
 * @date
 * @author PoKang Chen
 * @brief 
 */


#include <gtsam/sam/BearingRangeFactor.h>
#include <gtsam/slam/BetweenFactor.h>
#include <gtsam/slam/dataset.h>
#include <gtsam/geometry/Point3.h>
#include <gtsam/geometry/Pose2.h>
#include <gtsam/geometry/Rot3.h>
#include <gtsam/inference/FactorGraph.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/NonlinearFactor.h>
#include <gtsam/nonlinear/Values.h>
#include <gtsam/nonlinear/Values-inl.h>
#include <gtsam/nonlinear/GaussNewtonOptimizer.h>
#include <gtsam/linear/Sampler.h>
#include <gtsam/base/GenericValue.h>
#include <gtsam/base/Lie.h>
#include <gtsam/base/Matrix.h>
#include <gtsam/base/types.h>
#include <gtsam/base/Value.h>
#include <gtsam/base/Vector.h>
#include <gtsam/slam/PriorFactor.h>
#include <ReadNoiseModel.h>
#include <boost/assign/list_inserter.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <gtsam/nonlinear/ISAM2.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>
using namespace std;
using namespace gtsam;
namespace fs = boost::filesystem;

typedef std::pair<NonlinearFactorGraph::shared_ptr, Values::shared_ptr> GraphAndValues;
#define LINESIZE 81920
#define DATASET_DIR "@CMAKE_CURRENT_SOURCE_DIR@/dataset"
// #define GTSAM_INSTALLED_DATASET_DIR "@GTSAM_TOOLBOX_INSTALL_PATH@/gtsam_examples/Data"
/*
  Function for 
*/
boost::optional<IndexedPose> pVertex(istream& is, const string& tag) {
  if ((tag == "VERTEX2") || (tag == "VERTEX_SE2") || (tag == "VERTEX")) {
    Key id;
    double x, y, yaw;
    is >> id >> x >> y >> yaw;
    return IndexedPose(id, Pose2(x, y, yaw));
  } else {
    return boost::none;
  }
}
boost::optional<IndexedEdge> pEdge(istream& is, const string& tag) {
  if ((tag == "EDGE2") || (tag == "EDGE") || (tag == "EDGE_SE2")
      || (tag == "ODOMETRY")) {

    Key id1, id2;
    double x, y, yaw;
    is >> id1 >> id2 >> x >> y >> yaw;
    return IndexedEdge(pair<Key, Key>(id1, id2), Pose2(x, y, yaw));
  } else {
    return boost::none;
  }
}

std::pair<vector<boost::optional<IndexedPose>>, Values::shared_ptr> readVertex(){
    ifstream is;
    is.open("../dataset/input_INTEL_g2o.txt");
    if (!is)
    throw invalid_argument("load2D: can not find file input_INTEL_g2o");
    vector<boost::optional<IndexedPose>> vec;
    Values::shared_ptr initial(new Values);
    NonlinearFactorGraph::shared_ptr graph(new NonlinearFactorGraph);
    //------------------------
    string tag;
    int maxID = 0;
    bool addNoise = false;
    bool smart = 1;
    SharedNoiseModel model;
    //------------------------
    while (!is.eof()) {
        if (!(is >> tag))
        break;
        const auto indexed_pose = pVertex(is, tag);
        if (indexed_pose) {
        Key id = indexed_pose->first;

        // optional filter
        if (maxID && id >= maxID)
            continue;
        vec.push_back(indexed_pose);
        initial->insert(id, indexed_pose->second);

        }
        is.ignore(LINESIZE, '\n');
    }
    is.clear(); /* clears the end-of-file and error flags */
    is.seekg(0, ios::beg);
    return make_pair(vec, initial); 
    //vec: indexed_pose
    //Initial: id + indexed_pose
    
}
std::pair <vector<boost::optional<IndexedEdge>> ,vector<NonlinearFactor::shared_ptr>> readEdge(){
    ifstream is;
    is.open("../dataset/input_INTEL_g2o.txt");
    if (!is)
    throw invalid_argument("load2D: can not find file input_INTEL_g2o");
    vector<boost::optional<IndexedPose>> vec;
    Values::shared_ptr initial(new Values);
    NonlinearFactorGraph::shared_ptr graph(new NonlinearFactorGraph);
    //------------------------
    string tag;
    int maxID = 0;
    bool addNoise = false;
    bool smart = 1;
    SharedNoiseModel model;
    //------------------------
    while (!is.eof()) {
        if (!(is >> tag))
        break;

        const auto indexed_pose = pVertex(is, tag);
        if (indexed_pose) {
        Key id = indexed_pose->first;

        // optional filter
        if (maxID && id >= maxID)
            continue;
        
        // vec.push_back(indexed_pose);
        }
        is.ignore(LINESIZE, '\n');
    }
    is.clear(); /* clears the end-of-file and error flags */
    is.seekg(0, ios::beg);
    Sampler sampler;
    if (addNoise) {
        noiseModel::Diagonal::shared_ptr noise;
        if (model)
        noise = boost::dynamic_pointer_cast<noiseModel::Diagonal>(model);
        if (!noise)
        throw invalid_argument(
            "gtsam::load2D: invalid noise model for adding noise"
                "(current version assumes diagonal noise model)!");
        sampler = Sampler(noise);
    }
    Key id1, id2;
    bool haveLandmark = false;
    const bool useModelInFile = !model;
    vector<boost::optional<IndexedEdge>> vec2;
    vector<NonlinearFactor::shared_ptr> v_information;
    while (!is.eof()) {
        if (!(is >> tag))
        break;
        auto between_pose = pEdge(is, tag);
        if (between_pose) {
        vec2.push_back(between_pose); /// Added //////

        std::tie(id1, id2) = between_pose->first;
        Pose2& l1Xl2 = between_pose->second;

        // read noise model
        SharedNoiseModel modelInFile = readNoiseModel(is, smart, NoiseFormatG2O,
            KernelFunctionTypeNONE);

        // optional filter
        if (maxID && (id1 >= maxID || id2 >= maxID))
            continue;

        if (useModelInFile)
            model = modelInFile;

        if (addNoise)
            l1Xl2 = l1Xl2.retract(sampler.sample());

        // Insert vertices if pure odometry file
        if (!initial->exists(id1))
            initial->insert(id1, Pose2());
        if (!initial->exists(id2))
            initial->insert(id2, initial->at<Pose2>(id1) * l1Xl2);

        NonlinearFactor::shared_ptr factor(
            new BetweenFactor<Pose2>(id1, id2, l1Xl2, model));
        
        v_information.push_back(factor);
        graph->push_back(factor);//////Add graph

        }
        
    }

    return make_pair(vec2, v_information);
}


//  Writeg2o
// ----------------------------------------------------------------------------
void writeg2o(const NonlinearFactorGraph& graph, const Values& estimate,
    const string& filename) {
  fstream stream(filename.c_str(), fstream::out);

  // save 2D & 3D poses
  for (const auto& key_value : estimate) {
    auto p = dynamic_cast<const GenericValue<Pose2>*>(&key_value.value);
    if (!p) continue;
    const Pose2& pose = p->value();
    stream << "VERTEX_SE2 " << key_value.key << " " << pose.x() << " "
        << pose.y() << " " << pose.theta() << endl;
  }

  for(const auto& key_value: estimate) {
      auto p = dynamic_cast<const GenericValue<Pose3>*>(&key_value.value);
      if (!p) continue;
      const Pose3& pose = p->value();
      Point3 t = pose.translation();
      Rot3 R = pose.rotation();
      stream << "VERTEX_SE3:QUAT " << key_value.key << " " << t.x() << " "  << t.y() << " " << t.z()
        << " " << R.toQuaternion().x() << " " << R.toQuaternion().y() << " " << R.toQuaternion().z()
        << " " << R.toQuaternion().w() << endl;
  }

  // save edges (2D or 3D)
  for(const auto& factor_: graph) {
    boost::shared_ptr<BetweenFactor<Pose2> > factor =
        boost::dynamic_pointer_cast<BetweenFactor<Pose2> >(factor_);
    if (factor){
      SharedNoiseModel model = factor->noiseModel();
      boost::shared_ptr<noiseModel::Gaussian> gaussianModel =
          boost::dynamic_pointer_cast<noiseModel::Gaussian>(model);
      if (!gaussianModel){
        model->print("model\n");
        throw invalid_argument("writeG2o: invalid noise model!");
      }
      Matrix Info = gaussianModel->R().transpose() * gaussianModel->R();
      Pose2 pose = factor->measured(); //.inverse();
      stream << "EDGE_SE2 " << factor->key1() << " " << factor->key2() << " "
          << pose.x() << " " << pose.y() << " " << pose.theta();
      for (int i = 0; i < 3; i++){
        for (int j = i; j < 3; j++){
          stream << " " << Info(i, j);
        }
      }
      stream << endl;
    }

    boost::shared_ptr< BetweenFactor<Pose3> > factor3D =
        boost::dynamic_pointer_cast< BetweenFactor<Pose3> >(factor_);

    if (factor3D){
      SharedNoiseModel model = factor3D->noiseModel();

      boost::shared_ptr<noiseModel::Gaussian> gaussianModel =
          boost::dynamic_pointer_cast<noiseModel::Gaussian>(model);
      if (!gaussianModel){
        model->print("model\n");
        throw invalid_argument("writeG2o: invalid noise model!");
      }
      Matrix Info = gaussianModel->R().transpose() * gaussianModel->R();
      Pose3 pose3D = factor3D->measured();
      Point3 p = pose3D.translation();
      Rot3 R = pose3D.rotation();

      stream << "EDGE_SE3:QUAT " << factor3D->key1() << " " << factor3D->key2() << " "
          << p.x() << " "  << p.y() << " " << p.z()  << " " << R.toQuaternion().x()
          << " " << R.toQuaternion().y() << " " << R.toQuaternion().z()  << " " << R.toQuaternion().w();

      Matrix InfoG2o = I_6x6;
      InfoG2o.block(0,0,3,3) = Info.block(3,3,3,3); // cov translation
      InfoG2o.block(3,3,3,3) = Info.block(0,0,3,3); // cov rotation
      InfoG2o.block(0,3,3,3) = Info.block(0,3,3,3); // off diagonal
      InfoG2o.block(3,0,3,3) = Info.block(3,0,3,3); // off diagonal

      for (int i = 0; i < 6; i++){
        for (int j = i; j < 6; j++){
          stream << " " << InfoG2o(i, j);
        }
      }
      stream << endl;
    }
  }
  stream.close();
}


// 3D Case
// ----------------------------------------------------------------------------
// Parse3DFactors
BetweenFactorPose3s p3DFactors(const string& filename) {
  ifstream is(filename.c_str());
  if (!is) throw invalid_argument("parse3DFactors: can not find file " + filename);

  std::vector<BetweenFactor<Pose3>::shared_ptr> factors;
  while (!is.eof()) {
    char buf[LINESIZE];
    is.getline(buf, LINESIZE);
    istringstream ls(buf);
    string tag;
    ls >> tag;

    if (tag == "EDGE3") {
      Key id1, id2;
      double x, y, z, roll, pitch, yaw;
      ls >> id1 >> id2 >> x >> y >> z >> roll >> pitch >> yaw;
      Matrix m(6, 6);
      for (size_t i = 0; i < 6; i++)
        for (size_t j = i; j < 6; j++) ls >> m(i, j);
      SharedNoiseModel model = noiseModel::Gaussian::Information(m);
      factors.emplace_back(new BetweenFactor<Pose3>(
          id1, id2, Pose3(Rot3::Ypr(yaw, pitch, roll), {x, y, z}), model));
    }
    if (tag == "EDGE_SE3:QUAT") {
      Key id1, id2;
      double x, y, z, qx, qy, qz, qw;
      ls >> id1 >> id2 >> x >> y >> z >> qx >> qy >> qz >> qw;
      Matrix m(6, 6);
      for (size_t i = 0; i < 6; i++) {
        for (size_t j = i; j < 6; j++) {
          double mij;
          ls >> mij;
          m(i, j) = mij;
          m(j, i) = mij;
        }
      }
      Matrix mgtsam(6, 6);

      mgtsam.block<3, 3>(0, 0) = m.block<3, 3>(3, 3);  // cov rotation
      mgtsam.block<3, 3>(3, 3) = m.block<3, 3>(0, 0);  // cov translation
      mgtsam.block<3, 3>(0, 3) = m.block<3, 3>(0, 3);  // off diagonal
      mgtsam.block<3, 3>(3, 0) = m.block<3, 3>(3, 0);  // off diagonal

      SharedNoiseModel model = noiseModel::Gaussian::Information(mgtsam);
      factors.emplace_back(new BetweenFactor<Pose3>(
          id1, id2, Pose3(Rot3::Quaternion(qw, qx, qy, qz), {x, y, z}), model));
    }
  }
  return factors;
}

// Parse3Poses
std::map<Key, Pose3> p3DPoses(const string& filename) {
  ifstream is(filename.c_str());
  if (!is)
    throw invalid_argument("parse3DPoses: can not find file " + filename);

  std::map<Key, Pose3> poses;
  while (!is.eof()) {
    char buf[LINESIZE];
    is.getline(buf, LINESIZE);
    istringstream ls(buf);
    string tag;
    ls >> tag;

    if (tag == "VERTEX3") {
      Key id;
      double x, y, z, roll, pitch, yaw;
      ls >> id >> x >> y >> z >> roll >> pitch >> yaw;
      poses.emplace(id, Pose3(Rot3::Ypr(yaw, pitch, roll), {x, y, z}));
    }
    if (tag == "VERTEX_SE3:QUAT") {
      Key id;
      double x, y, z, qx, qy, qz, qw;
      ls >> id >> x >> y >> z >> qx >> qy >> qz >> qw;
      poses.emplace(id, Pose3(Rot3::Quaternion(qw, qx, qy, qz), {x, y, z}));
    }
  }
  return poses;
}


// ------------------------- For Iteration Method
// Parse3DFactors
BetweenFactorPose3s p3DFactors_3D_increments(const string& filename) {
  ifstream is(filename.c_str());
  if (!is) throw invalid_argument("parse3DFactors: can not find file " + filename);

  std::vector<BetweenFactor<Pose3>::shared_ptr> factors;
  while (!is.eof()) {
    char buf[LINESIZE];
    is.getline(buf, LINESIZE);
    istringstream ls(buf);
    string tag;
    ls >> tag;

    if (tag == "EDGE3") {
      Key id1, id2;
      double x, y, z, roll, pitch, yaw;
      ls >> id1 >> id2 >> x >> y >> z >> roll >> pitch >> yaw;
      Matrix m(6, 6);
      for (size_t i = 0; i < 6; i++)
        for (size_t j = i; j < 6; j++) ls >> m(i, j);
      SharedNoiseModel model = noiseModel::Gaussian::Information(m);
      factors.emplace_back(new BetweenFactor<Pose3>(
          id1, id2, Pose3(Rot3::Ypr(yaw, pitch, roll), {x, y, z}), model));
    }
    if (tag == "EDGE_SE3:QUAT") {
      Key id1, id2;
      double x, y, z, qx, qy, qz, qw;
      ls >> id1 >> id2 >> x >> y >> z >> qx >> qy >> qz >> qw;
      Matrix m(6, 6);
      for (size_t i = 0; i < 6; i++) {
        for (size_t j = i; j < 6; j++) {
          double mij;
          ls >> mij;
          m(i, j) = mij;
          m(j, i) = mij;
        }
      }
      Matrix mgtsam(6, 6);

      mgtsam.block<3, 3>(0, 0) = m.block<3, 3>(3, 3);  // cov rotation
      mgtsam.block<3, 3>(3, 3) = m.block<3, 3>(0, 0);  // cov translation
      mgtsam.block<3, 3>(0, 3) = m.block<3, 3>(0, 3);  // off diagonal
      mgtsam.block<3, 3>(3, 0) = m.block<3, 3>(3, 0);  // off diagonal

      SharedNoiseModel model = noiseModel::Gaussian::Information(mgtsam);
      factors.emplace_back(new BetweenFactor<Pose3>(
          id1, id2, Pose3(Rot3::Quaternion(qw, qx, qy, qz), {x, y, z}), model));
      
    }
  }
  return factors;
}

