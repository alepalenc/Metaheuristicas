#include<iostream>
#include<vector>
#include<unordered_set>
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




double evaluateFitness(const vector<unsigned> &solution, const vector<vector<double> > &matrix){
    double fitness = 0.0;
    for(unsigned i=0; i<solution.size(); ++i){
        for(unsigned j=i+1; j<solution.size(); ++j)
            fitness += matrix[solution[i]][solution[j]];
    }
    return fitness;
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




void greedy(const vector<vector<double> > &matrix, unsigned m){
    clock_t start_time, elapsed;
    start_time = clock();
    
    vector<unsigned> sel;
    unordered_set<unsigned> not_sel;
    for(unsigned i=0; i<matrix.size(); ++i)
        not_sel.insert(i);
    
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
    
    sel.push_back(*it_best);
    not_sel.erase(it_best);
    
    while(sel.size() < m){
        best_dist = 0.0;
        for(it=not_sel.begin(); it!=not_sel.end(); ++it){
            dist = distToSel(matrix, *it, sel);
            if(dist > best_dist){
                best_dist = dist;
                it_best = it;
            }
        }
    
        sel.push_back(*it_best);
        not_sel.erase(it_best);
    }
    
    elapsed = clock() - start_time;
    
    cout << fixed;
    cout << evaluateFitness(sel, matrix) << "\t" << (double) elapsed / CLOCKS_PER_SEC << endl;
}




int main(int argc, char **argv){
    unsigned n, m;
    
    cin >> n >> m;
    vector<double> v(n, 0);
    vector<vector<double> > matrix(n, v);
    readInput(matrix);
    
    greedy(matrix, m);
    
    return 0;
}
