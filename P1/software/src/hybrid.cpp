#include<iostream>
#include<vector>
#include<set>
#include<unordered_set>
#include<utility> // pair
#include<algorithm> // random_shuffle
#include<cstdlib> // rand, srand
#include<ctime> // clock


using namespace std;


void readInput(vector<vector<double> > &matrix){
    double dist;
    for(unsigned i=0; i<matrix.size(); ++i){
        for(unsigned j=i+1; j<matrix.size(); ++j){
            cin >> dist >> dist >> dist;
            matrix[i][j] = matrix[j][i] = dist;
        }
    }
    for(unsigned i=0; i<matrix.size(); ++i)
        matrix[i][i] = 0.0;
}




struct contribution_comp{
    bool operator()(const pair<unsigned,double> &p1, const pair<unsigned,double> &p2) const{
        if(p1.second != p2.second)
            return p1.second < p2.second;
        else
            return p1.first < p2.first;
    }
};




double evaluateFitness(const set<pair<unsigned, double>, contribution_comp> &solution){
    double fitness = 0.0;
    set<pair<unsigned, double>, contribution_comp>::iterator it;
    for(it=solution.begin(); it!=solution.end(); ++it)
        fitness += it->second;
    return fitness/2;
}




double acumDist(const vector<vector<double> > &matrix, unsigned element){
    double acum_dist = 0.0;
    for(unsigned i=0; i<matrix.size(); ++i)
        acum_dist += matrix[element][i];
    return acum_dist;
}




double distToSel(const vector<vector<double> > &matrix, unsigned element, const vector<unsigned> &sel){
    double min_dist, dist;
    min_dist = matrix[element][sel[0]];
    for(unsigned i=1; i<sel.size(); ++i){
        dist = matrix[element][sel[i]];
        if(dist < min_dist)
            min_dist = dist;
    }
    return min_dist;
}




void greedy(set<pair<unsigned, double>, contribution_comp> &sel, 
            unordered_set<unsigned> &not_sel, 
            const vector<vector<double> > &matrix, 
            unsigned m){
    
    sel.clear();
    not_sel.clear();
    
    for(unsigned i=0; i<matrix.size(); ++i)
        not_sel.insert(i);
    
    vector<unsigned> sel_aux;
    unordered_set<unsigned>::iterator it_best, it;
    double dist, best_dist;
    best_dist = 0.0;
    for(it=not_sel.begin(); it!=not_sel.end(); ++it){
        dist = acumDist(matrix, *it);
        if(dist > best_dist){
            best_dist = dist;
            it_best = it;
        }
    }
    
    sel_aux.push_back(*it_best);
    not_sel.erase(it_best);
    
    while(sel_aux.size() < m){
        best_dist = 0.0;
        for(it=not_sel.begin(); it!=not_sel.end(); ++it){
            dist = distToSel(matrix, *it, sel_aux);
            if(dist > best_dist){
                best_dist = dist;
                it_best = it;
            }
        }
    
        sel_aux.push_back(*it_best);
        not_sel.erase(it_best);
    }
    
    vector<double> contribs(m, 0.0);
    for(unsigned i=0; i<sel_aux.size(); ++i){
        for(unsigned j=i+1; j<sel_aux.size(); ++j){
            contribs[i] += matrix[sel_aux[i]][sel_aux[j]];
            contribs[j] += matrix[sel_aux[i]][sel_aux[j]];
        }
        sel.insert(make_pair(sel_aux[i], contribs[i]));
    }
}





double calculateContribution(unsigned element, 
                             unsigned element_to_ignore, 
                             const set<pair<unsigned, double>, contribution_comp> &sel, 
                             const vector<vector<double> > &matrix){
    
    double contribution = -matrix[element][element_to_ignore];
    set<pair<unsigned, double>, contribution_comp>::iterator it;
    
    for(it=sel.begin(); it!=sel.end(); ++it)
        contribution += matrix[element][it->first];
    
    return contribution;
}





void exchangeNeighbour(set<pair<unsigned, double>, contribution_comp> &sel, 
                       unordered_set<unsigned> &not_sel, 
                       const vector<vector<double> > &matrix, 
                       set<pair<unsigned, double>, contribution_comp>::iterator it_old_element, 
                       unordered_set<unsigned>::iterator it_new_element, 
                       double new_element_contrib){
    
    set<pair<unsigned, double>, contribution_comp> sel_aux;
    set<pair<unsigned, double>, contribution_comp>::iterator it_sel;
    unsigned old_element;
    double updated_contrib;
    
    old_element = it_old_element->first;
    sel.erase(it_old_element);
    not_sel.insert(old_element);
    
    for(it_sel=sel.begin(); it_sel!=sel.end(); ++it_sel){
        updated_contrib = it_sel->second + matrix[it_sel->first][*it_new_element] - matrix[it_sel->first][old_element];
        sel_aux.insert(make_pair(it_sel->first, updated_contrib));
    }
    
    sel.swap(sel_aux);
    
    sel.insert(make_pair(*it_new_element, new_element_contrib));
    not_sel.erase(it_new_element);
}






void localSearch(const vector<vector<double> > &matrix, unsigned m){
    clock_t start_time, elapsed;
    start_time = clock();
    
    set<pair<unsigned, double>, contribution_comp> sel;
    unordered_set<unsigned> not_sel;
    set<pair<unsigned, double>, contribution_comp>::iterator it_sel;
    unordered_set<unsigned>::iterator it_not_sel;
    double new_element_contrib;
    bool exchange;
    
    greedy(sel, not_sel, matrix, m);
    
    do{
        exchange = false;
        for(it_sel=sel.begin(); it_sel!=sel.end() and !exchange; ++it_sel){
            for(it_not_sel=not_sel.begin(); it_not_sel!=not_sel.end() and !exchange; ++it_not_sel){
                new_element_contrib = calculateContribution(*it_not_sel, it_sel->first, sel, matrix);
                if(new_element_contrib > it_sel->second){
                    exchangeNeighbour(sel, not_sel, matrix, it_sel, it_not_sel, new_element_contrib);
                    exchange = true;
                }
            }
        }
    }while(exchange);
    
    elapsed = clock() - start_time;
    
    cout << fixed;
    cout << evaluateFitness(sel) << "\t" << (double) elapsed / CLOCKS_PER_SEC << endl;
}




int main(int argc, char **argv){
    unsigned n, m;
    
    cin >> n >> m;
    vector<double> v(n, 0);
    vector<vector<double> > matrix(n, v);
    readInput(matrix);
    
    srand(1);
    localSearch(matrix, m);
    
    return 0;
}
