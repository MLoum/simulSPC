#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <sstream>

struct SpeciesConfig {
    double concentration_nM{5.0};
    double r_hydro_nm{5.0};
    double r_hydro_dispersion_pct{0.0};
    double lifetime_ns{3.8};
    double lifetime_dispersion_pct{0.0};
    double cross_section_cm2{3.0e-16};
    double quantum_yield{0.8};
    bool is_isotropic{false};
};

class Config {
private:
    std::unordered_map<std::string, std::string> globals_;
    std::vector<SpeciesConfig> species_list_;

    static inline std::string trim(const std::string& s) {
        size_t start = s.find_first_not_of(" \t\r\n");
        if (start == std::string::npos) return "";
        size_t end = s.find_last_not_of(" \t\r\n");
        return s.substr(start, end - start + 1);
    }

public:
    Config() = default;

    bool load(const std::string& filepath) {
        std::ifstream file(filepath);
        if (!file.is_open()) return false;

        std::string line;
        bool in_species_block = false;
        SpeciesConfig current_species;

        while (std::getline(file, line)) {
            size_t comment_pos = line.find('#');
            if (comment_pos != std::string::npos) {
                line = line.substr(0, comment_pos);
            }
            line = trim(line);
            if (line.empty()) continue;

            if (line == "[Species]") {
                if (in_species_block) {
                    species_list_.push_back(current_species);
                }
                in_species_block = true;
                current_species = SpeciesConfig{};
                continue;
            }

            size_t sep = line.find('=');
            if (sep != std::string::npos) {
                std::string key = trim(line.substr(0, sep));
                std::string val = trim(line.substr(sep + 1));

                if (in_species_block) {
                    if (key == "concentration_nM") current_species.concentration_nM = std::stod(val);
                    else if (key == "r_hydro_nm") current_species.r_hydro_nm = std::stod(val);
                    else if (key == "r_hydro_dispersion_pct") current_species.r_hydro_dispersion_pct = std::stod(val);
                    else if (key == "lifetime_ns") current_species.lifetime_ns = std::stod(val);
                    else if (key == "lifetime_dispersion_pct") current_species.lifetime_dispersion_pct = std::stod(val);
                    else if (key == "cross_section_cm2") current_species.cross_section_cm2 = std::stod(val);
                    else if (key == "quantum_yield") current_species.quantum_yield = std::stod(val);
                    else if (key == "is_isotropic") current_species.is_isotropic = (val == "1" || val == "true");
                } else {
                    globals_[key] = val;
                }
            }
        }

        if (in_species_block) {
            species_list_.push_back(current_species);
        }
        return true;
    }

    double get_double(const std::string& key, double default_val) const {
        auto it = globals_.find(key);
        if (it != globals_.end()) {
            try { return std::stod(it->second); } catch (...) {}
        }
        return default_val;
    }

    int get_int(const std::string& key, int default_val) const {
        auto it = globals_.find(key);
        if (it != globals_.end()) {
            try { return std::stoi(it->second); } catch (...) {}
        }
        return default_val;
    }

    const std::vector<SpeciesConfig>& species() const { return species_list_; }
};