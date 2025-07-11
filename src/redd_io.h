#ifndef REDD_IO_H
#define REDD_IO_H

#include <highfive/H5File.hpp>

#include "f5c.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdio>
#include <unordered_map>

inline void init_redd_candidate_ratio_map(core_t* core){
    core->redd_candidate_ratio_map = std::unordered_map<std::string, std::unordered_map<u_int64_t, float>>();
    std::ifstream file(core->opt.redd_candidate_file);
    if (!file.is_open()) {
        std::cerr << "Error opening file\n";
    }

    std::string line;
    while (std::getline(file, line)) {
        std::vector<std::string> columns;
        std::istringstream iss(line);
        std::string cell;

        while (std::getline(iss, cell, '\t')) {
            columns.push_back(cell);
        }
        std::string contig = columns[0];
        // candidate file is 1-based index, use 0-based index in program
        u_int64_t ref_position = std::stoi(columns[1]) - 1;
        float ratio = std::stof(columns[3]);
        if (core->redd_candidate_ratio_map.find(contig) == core->redd_candidate_ratio_map.end()){
            core->redd_candidate_ratio_map[contig] = std::unordered_map<u_int64_t, float>();
        }
        core->redd_candidate_ratio_map[contig][ref_position] = ratio;

    }

    file.close();


}
inline void init_redd_hdf5_file(core_t* core,std::string output_file){
    std::remove(output_file.c_str());
    HighFive::File file(output_file, HighFive::File::ReadWrite | HighFive::File::Create);
    HighFive::DataSetCreateProps X_props;
    X_props.add(HighFive::Chunking({10240,core->opt.redd_window_size,5}));
    HighFive::DataSetCreateProps y_props;
    y_props.add(HighFive::Chunking({10240,core->opt.redd_window_size}));
    HighFive::DataSetCreateProps ratio_props;
    ratio_props.add(HighFive::Chunking({10240}));
    // Create a dataset with an unlimited dimension for appending
    file.createDataSet<float>(
        "/X",
        HighFive::DataSpace({0,core->opt.redd_window_size,5}, {HighFive::DataSpace::UNLIMITED,core->opt.redd_window_size,5}), // Initial size 10, unlimited max
        X_props // Chunking required for extendable datasets
    );
    file.createDataSet<uint32_t>(
        "/y_ref",
        HighFive::DataSpace({0,core->opt.redd_window_size}, {HighFive::DataSpace::UNLIMITED,core->opt.redd_window_size}), // Initial size 10, unlimited max
        y_props // Chunking required for extendable datasets
    );
    file.createDataSet<uint32_t>(
        "/y_call",
        HighFive::DataSpace({0,core->opt.redd_window_size}, {HighFive::DataSpace::UNLIMITED,core->opt.redd_window_size}), // Initial size 10, unlimited max
        y_props // Chunking required for extendable datasets
    );
    file.createDataSet<float>(
        "/ratio",
        HighFive::DataSpace({0}, {HighFive::DataSpace::UNLIMITED}), // Initial size 10, unlimited max
        ratio_props // Chunking required for extendable datasets
    );
};
inline void append_arr_to_dataset(core_t* core,std::string output_file,
    std::vector<std::vector<std::vector<float>>> X_arr, 
    std::vector<std::vector<uint32_t>> y_ref_arr,
    std::vector<std::vector<uint32_t>> y_call_arr,
    std::vector<float> ratio_arr){
    HighFive::File hdf5_output(output_file, HighFive::File::ReadWrite);
    auto ratio_dataset = hdf5_output.getDataSet("ratio");
    auto X_dataset = hdf5_output.getDataSet("X");
    auto y_ref_dataset = hdf5_output.getDataSet("y_ref");
    auto y_call_dataset = hdf5_output.getDataSet("y_call");
    size_t existing_shape = X_dataset.getDimensions()[0];
    X_dataset.resize({existing_shape + X_arr.size(), core->opt.redd_window_size,5});
    X_dataset.select({existing_shape, 0,0}, {X_arr.size(), core->opt.redd_window_size, 5}).write(X_arr);
    existing_shape = y_ref_dataset.getDimensions()[0];
    y_ref_dataset.resize({ existing_shape + y_ref_arr.size(), core->opt.redd_window_size});
    y_ref_dataset.select({existing_shape, 0}, {y_ref_arr.size(), core->opt.redd_window_size}).write(y_ref_arr);
    existing_shape = y_call_dataset.getDimensions()[0];
    y_call_dataset.resize({existing_shape + y_call_arr.size(), core->opt.redd_window_size});
    y_call_dataset.select({existing_shape, 0}, {y_call_arr.size(), core->opt.redd_window_size}).write(y_call_arr);
    if (!ratio_arr.empty()){
        existing_shape = ratio_dataset.getDimensions()[0];
        ratio_dataset.resize({existing_shape + ratio_arr.size()});
        ratio_dataset.select({existing_shape}, {ratio_arr.size()}).write(ratio_arr);

    }


};
#endif
