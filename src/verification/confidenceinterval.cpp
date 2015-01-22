#include "confidenceinterval.hpp"

void VerifyAttack::confidenceinterval::init(){
    if (alphaArg.getValue() <= 0 || alphaArg.getValue() >= 1){
        cerr << "Alpha must be 0 < alpha < 1" << endl;
        exit ( 1 );
    }
}

void VerifyAttack::confidenceinterval::batchBest( unsigned long long id, shared_ptr<StatisticIndexMatrix>& s ) {
    Eigen::Matrix < StatisticValueType, 1, KEYNUM > maxPearson;
    ( *s ) = s->cwiseAbs();
    maxPearson = s->colwise().maxCoeff();
    checkMutex.lock();
       for(unsigned long long key = 0; key < KEYNUM; key++){
           //TODO: Maybe transform the two vectors in a 2xKEYNUM Matrix, compute the colwise max and assign it to bestPearson
           if (bestPearson (key) < maxPearson (key)){
               bestPearson (key) = maxPearson (key);
           }
       }
    checkMutex.unlock();
}

bool VerifyAttack::confidenceinterval::verify(shared_ptr<StatisticProg::base> stat){
    unsigned long long maxRow, maxCol;
    StatisticValueType best = bestPearson.maxCoeff(&maxRow, &maxCol);
    bestPearson (maxCol) = 0;
    ConfidencePair bestConf = getConfidence(best, currentTraces, alphaArg.getValue());

    StatisticValueType secondBest = bestPearson.maxCoeff();
    ConfidencePair secondBestConf = getConfidence(secondBest, currentTraces, alphaArg.getValue());
    
    bestPearson = Trace::Zero(KEYNUM);
    if (bestConf.first > secondBestConf.second){
        return true;
    } else {
        return false;
    }
}


ConfidencePair VerifyAttack::confidenceinterval::getConfidence(StatisticValueType r, unsigned long n, StatisticValueType alpha)
{
    //TODO: if n < 4 I think there will be problems

    // Confidence level
    StatisticValueType gamma = 1-alpha;
    // Fisher transform r to map correctly the domains
    StatisticValueType Fr = atanh(r);
    StatisticValueType z_score = sqrt(n-3)*Fr;
    // obtain Gaussian quantile of transformed z_score
    StatisticValueType half_width=sqrt(2)*boost::math::erf_inv(2*((1+gamma)/2)-1);
    // compute confidence interval in Fisher domain
    StatisticValueType inf=z_score-half_width;
    StatisticValueType sup=z_score+half_width;
    // invert Fisher map to get the domain straight again
    inf=tanh(inf/sqrt(n-3));
    sup=tanh(sup/sqrt(n-3));
    return ConfidencePair(inf,sup);
}
