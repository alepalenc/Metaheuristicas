#include<iostream>
#include<vector>
#include<cmath> // log, exp
#include<algorithm> // swap, random_shuffle
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




struct Solution{
    vector<unsigned> selected;
    vector<unsigned> not_selected;
    double fitness = -1.0;
};




void evaluateFitness(Solution &sol, const vector<vector<double> > &matrix){
    sol.fitness = 0.0;
    
    for(unsigned i=0; i<sol.selected.size(); ++i){
        for(unsigned j=i+1; j<sol.selected.size(); ++j)
            sol.fitness += matrix[sol.selected[i]][sol.selected[j]];
    }
}




void generateRandomSolution(Solution &sol, unsigned n, unsigned m){
    unsigned random_element;
    
    sol = Solution();
    
    for(unsigned i=0; i<n; ++i)
        sol.not_selected.push_back(i);
    
    random_shuffle(sol.not_selected.begin(), sol.not_selected.end());
    
    for(unsigned i=0; i<m; ++i){
        // genero entero entre 0 y sol.not_selected.size()-1
        random_element = rand() % sol.not_selected.size();
        // intercambio los valores de sol.not_selected.back() y sol.not_selected[random_element]
        swap(sol.not_selected.back(), sol.not_selected[random_element]);
        // añado como gen el último elemento de not_selected
        sol.selected.push_back(sol.not_selected.back());
        // elimino el último elemento de not_selected
        sol.not_selected.pop_back();
    }
}




double contribution(const Solution &sol, unsigned element, const vector<vector<double> > &matrix){
    double contrib = 0.0;
    
    for(unsigned i=0; i<sol.selected.size(); ++i)
        contrib += matrix[element][sol.selected[i]];
    
    return contrib;
}





void ES(const vector<vector<double> > &matrix, unsigned m){
    clock_t start_time, elapsed;
    start_time = clock();
    
    Solution sol, best_sol;
    generateRandomSolution(sol, matrix.size(), m);
    evaluateFitness(sol, matrix);
    best_sol = sol;
    
    const unsigned MAX_VECINOS = 10*m,
                   MAX_EXITOS = m,
                   MAX_EVALUATIONS = 100000,
                   M = MAX_EVALUATIONS/MAX_VECINOS;
    const double PHI = 0.3,
                 MU = 0.3,
                 T0 = (MU*sol.fitness)/(-log(PHI));
    
    double t = 0.001;
    while(t >= T0)
        t *= 0.1;
    
    const double TF = t,
                 BETA = (T0-TF)/(M*T0*TF);
    
    t = T0;
    unsigned evaluations = 1;
    unsigned num_vecinos, num_exitos;
    unsigned sel_element, not_sel_element;
    double inc_fitness;
    
    do{
        num_vecinos = 0;
        num_exitos = 0;
        
        while( num_vecinos < MAX_VECINOS and num_exitos < MAX_EXITOS and evaluations < MAX_EVALUATIONS ){
            // elijo aleatoriamente un elemento seleccionado
            sel_element = rand() % sol.selected.size();
            // elijo aleatoriamente un elemento no seleccionado
            not_sel_element = rand() % sol.not_selected.size();
            // calculo el incremento del fitness entre la nueva solución y la actual
            inc_fitness = ( contribution(sol,sol.not_selected[not_sel_element],matrix)
                          - matrix[sol.selected[sel_element]][sol.not_selected[not_sel_element]] )
                          - contribution(sol,sol.selected[sel_element],matrix);
            
            if( inc_fitness > 0 or (1.0*rand())/RAND_MAX <= exp(inc_fitness/t) ){
                swap(sol.selected[sel_element], sol.not_selected[not_sel_element]);
                sol.fitness += inc_fitness;
                ++num_exitos;
                if( sol.fitness > best_sol.fitness )
                    best_sol = sol;
            }
            
            ++num_vecinos;
            ++evaluations;
        }
        
        t = t/(1+BETA*t);
    }while( evaluations < MAX_EVALUATIONS and num_exitos > 0 );
    
    
    elapsed = clock() - start_time;
    
    cout << fixed;
    cout << best_sol.fitness << "\t" << (double) elapsed / CLOCKS_PER_SEC << endl;
    
}




int main(int argc, char **argv){
    unsigned n, m;
    
    cin >> n >> m;
    vector<double> v(n, 0);
    vector<vector<double> > matrix(n, v);
    readInput(matrix);
    
    srand(1);
    ES(matrix, m);
    
    return 0;
}
