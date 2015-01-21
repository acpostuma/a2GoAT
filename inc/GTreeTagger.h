#ifndef __GTreeTagger_h__
#define __GTreeTagger_h__


#include <TLorentzVector.h>
#include <TClonesArray.h>


#include "GTree.h"


#define GTreeTagger_MAX 4096

class  GTreeTagger : public GTree
{
private:
    Int_t           nTagged;
    Int_t           taggedChannel[GTreeTagger_MAX];
    Double_t        taggedTime[GTreeTagger_MAX];
    Double_t        taggedEnergy[GTreeTagger_MAX];

protected:
    virtual void    SetBranchAdresses();
    virtual void    SetBranches();

public:
    GTreeTagger(GTreeManager *Manager);
    virtual ~GTreeTagger();

    virtual void    Clear()             {nTagged = 0;}

            Int_t           GetNTagged()                        const	{return nTagged;}
    const	Double_t*       GetTaggedEnergy()                   const	{return taggedEnergy;}
            Double_t        GetTaggedEnergy(const Int_t index)	const	{return taggedEnergy[index];}
    const	Int_t*          GetTaggedChannel()                  const	{return taggedChannel;}
            Int_t           GetTaggedChannel(const Int_t index) const	{return taggedChannel[index];}
    const	Double_t*       GetTaggedTime()                     const	{return taggedTime;}
            Double_t        GetTaggedTime(const Int_t index)    const	{return taggedTime[index];}
    TLorentzVector          GetVector(const Int_t index)        const   {return TLorentzVector(0, 0, taggedEnergy[index], taggedEnergy[index]);}
    TLorentzVector          GetVectorProtonTarget(const Int_t index)    const;
};

#endif
