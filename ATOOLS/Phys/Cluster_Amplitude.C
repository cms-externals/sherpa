#include "ATOOLS/Phys/Cluster_Amplitude.H"

#include <algorithm>
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Phys/Flow.H"
#include "ATOOLS/Org/Smart_Pointer.C"
#include "ATOOLS/Phys/Blob.H"


using namespace ATOOLS;

ClusterAmplitude_PVector::ClusterAmplitude_PVector()
{
}

ClusterAmplitude_PVector::~ClusterAmplitude_PVector()
{
  while (!empty()) {
    Cluster_Amplitude *ampl(back());
    pop_back();
    delete ampl;
  }
}

Cluster_Amplitude * ClusterAmplitude_PVector::OneAmpl()
{
  // links all amplitudes from this amplitude vector and return the linked, single amplitude
  if (this->size()==0) return Cluster_Amplitude::New();
  Cluster_Amplitude * tmp = this->at(0)->CopyAll();
  for(int i=1; i< this->size(); i++){
      if (this->at(i)){
          Cluster_Amplitude * next_prev = this->at(i)->CopyAll();
          while (next_prev->Prev()) next_prev = next_prev->Prev();
          tmp->SetNext(next_prev);
          while (tmp->Next()) tmp = tmp->Next();

        }
    }

  while (tmp->Prev()) tmp = tmp->Prev();
  return tmp;
}

ClusterAmplitude_PVector Cluster_Amplitude::s_ampls;

Cluster_Amplitude::Cluster_Amplitude(Cluster_Amplitude *const prev):
  p_prev(prev), p_next(NULL), 
  m_oew(0), m_oqcd(0), m_nin(0), m_new(0), m_ncl(0),
  m_kin(0), m_nlo(0), m_flag(0),
  m_mur2(0.0), m_muf2(0.0), m_muq2(0.0), m_mu2(0.0),
  m_kt2(0.0), m_z(0.0), m_phi(0.0), m_lkf(0.0),
  p_jf(NULL), p_proc(NULL), p_procs(NULL), p_dinfo(NULL), p_ms(NULL)
{
  if (p_prev!=NULL) p_prev->p_next=this;
}

Cluster_Amplitude::~Cluster_Amplitude()
{
  if (p_next) p_next->Delete();
  for (size_t i(0);i<m_legs.size();++i) m_legs[i]->Delete();
  if (p_prev) p_prev->p_next=NULL;
}

Cluster_Amplitude *Cluster_Amplitude::New
(Cluster_Amplitude *const prev)
{
  if (s_ampls.empty()) {
    return new Cluster_Amplitude(prev);
  }
  Cluster_Amplitude *ca(s_ampls.back());
  s_ampls.pop_back();
  ca->p_prev=prev;
  ca->p_next=NULL;
  ca->m_oew=ca->m_oqcd=0;
  ca->m_nin=ca->m_new=ca->m_ncl=ca->m_kin=ca->m_nlo=ca->m_flag=0;
  ca->m_mur2=ca->m_muf2=ca->m_muq2=ca->m_mu2=0.0;
  ca->m_kt2=ca->m_z=ca->m_phi=ca->m_lkf=0.0;
  ca->p_jf=ca->p_proc=ca->p_procs=ca->p_dinfo=NULL;
  ca->p_ms=NULL;
  if (ca->p_prev!=NULL) ca->p_prev->p_next=ca;
  return ca;
}

void Cluster_Amplitude::Delete()
{
  if (p_next) p_next->Delete();
  for (size_t i(0);i<m_legs.size();++i) m_legs[i]->Delete();
  m_legs.clear();
  m_decs.clear();
  m_cmap.clear();
  if (p_prev) p_prev->p_next=NULL;
  p_prev=p_next=NULL;
  s_ampls.push_back(this);
}

void Cluster_Amplitude::CreateLeg
(const Vec4D &p,const Flavour &fl,
 const ColorID &col,const size_t &id)
{
  m_legs.push_back(Cluster_Leg::New(this,p,fl,col));
  if (id!=std::string::npos) m_legs.back()->SetId(id);
  else m_legs.back()->SetId(1<<(m_legs.size()-1));
}

Cluster_Amplitude *Cluster_Amplitude::Copy() const
{
  Cluster_Amplitude *copy(Cluster_Amplitude::New());
  copy->CopyFrom(this);
  return copy;
}

void Cluster_Amplitude::CopyFrom
(const Cluster_Amplitude *const master,const int mode)
{
  Cluster_Amplitude *sprev(p_prev), *snext(p_next);
  *this=*master;
  p_prev=sprev;
  p_next=snext;
  if (mode==1) m_legs.clear();
  else for (size_t i(0);i<m_legs.size();++i)
    m_legs[i] = Cluster_Leg::New(this,*master->m_legs[i]);
}

Cluster_Amplitude *Cluster_Amplitude::CopyNext() const
{
  const Cluster_Amplitude *root(this);
  Cluster_Amplitude *prev(NULL), *ref(NULL);
  while (root) {
    Cluster_Amplitude *copy(root->Copy());
    if (prev!=NULL) (prev->p_next=copy)->p_prev=prev;
    prev=copy;
    if (root==this) ref=prev;
    root=root->Next();
  }
  return ref;
}

Cluster_Amplitude *Cluster_Amplitude::CopyAll() const
{
  const Cluster_Amplitude *root(this);
  Cluster_Amplitude *prev(NULL), *ref(NULL);
  while (root->Prev()) root=root->Prev();
  while (root) {
    Cluster_Amplitude *copy(root->Copy());
    if (prev!=NULL) (prev->p_next=copy)->p_prev=prev;
    prev=copy;
    if (root==this) ref=prev;
    root=root->Next();
  }
  return ref;
}

void Cluster_Amplitude::CombineLegs
(Cluster_Leg *const i,Cluster_Leg *const j,
 const Flavour &fl,const ColorID &col)
{
  if (i->Amplitude()!=this || j->Amplitude()!=this) 
    THROW(fatal_error,"Leg not owned by current amplitude");
  for (ClusterLeg_Vector::iterator clit(m_legs.begin());
       clit!=m_legs.end();++clit) {
    if (*clit==i || *clit==j) {
      *clit = Cluster_Leg::New(this,i->Mom()+j->Mom(),fl,col);
      (*clit)->SetId( (i->Id() | j->Id()) );
      i->Delete();
      j->Delete();
      for (++clit;clit!=m_legs.end();++clit)
	if (*clit==i || *clit==j) {
	  clit=m_legs.erase(clit);
	  break;
	}
      break;
    }
  }
}

Cluster_Amplitude *Cluster_Amplitude::InitNext()
{
  if (p_next!=NULL) p_next->Delete();
  p_next = Cluster_Amplitude::New(this);
  return p_next;
}

Cluster_Amplitude *Cluster_Amplitude::InitPrev()
{
  if (p_prev!=NULL) return NULL;
  p_prev = Cluster_Amplitude::New();
  p_prev->p_next=this;
  return p_prev;
}

void Cluster_Amplitude::SetNext(Cluster_Amplitude *const next) 
{
  if (p_next!=NULL) p_next->Delete();
  if (next->p_prev) next->p_prev->p_next=NULL;
  (p_next=next)->p_prev=this;
}

void Cluster_Amplitude::UnsetNext()
{
  if (p_next) p_next->p_prev=NULL;
  p_next=NULL;
}

void Cluster_Amplitude::DeletePrev()
{
  if (p_prev==NULL) return;
  p_prev->p_next=NULL;
  while (p_prev->Prev()) p_prev=p_prev->Prev();
  p_prev->Delete();
  p_prev=NULL;
}

void Cluster_Amplitude::DeleteNext()
{
  if (p_next!=NULL) p_next->Delete();
  p_next=NULL;
}

class Order_LegId {
public:
  int operator()(const Cluster_Leg *a,const Cluster_Leg *b) 
  { return a->Id()<b->Id(); }
};// end of class Order_LegId

void Cluster_Amplitude::IdSort()
{
  std::stable_sort(m_legs.begin(),m_legs.end(),Order_LegId());
}

size_t Cluster_Amplitude::NQCD() const
{
  size_t nqcd(0);
  for (size_t i(0);i<m_legs.size();++i)
    nqcd+=m_legs[i]->Flav().Strong();
  return nqcd;
}

size_t Cluster_Amplitude::NEW() const
{
  size_t nw(0);
  for (size_t i(0);i<m_legs.size();++i)
    nw+=!m_legs[i]->Flav().Strong();
  return nw;
}

Cluster_Leg *Cluster_Amplitude::IdLeg(const size_t &id) const
{
  for (size_t i(0);i<m_legs.size();i++)
    if (m_legs[i]->Id()==id) return m_legs[i];
  return NULL;
}

Cluster_Leg *Cluster_Amplitude::Splitter() const
{
  for (size_t i(0);i<m_legs.size();i++)
    if (m_legs[i]->K()) return m_legs[i];
  return NULL;
}

void Cluster_Amplitude::SetColours
(Cluster_Leg *const lij,Cluster_Leg *const li,Cluster_Leg *const lj)
{
  ColorID colij(lij->Col()), coli(0,0), colj(0,0);
  if (lij->Flav().StrongCharge()==3) {
    if (li->Flav().StrongCharge()==3) {
      if (lj->Flav().Strong()) {
        // triplet -> triplet octet
	size_t nc(Flow::Counter());
	colj.m_j=coli.m_i=nc;
	colj.m_i=colij.m_i;
      }
      else {
        // triplet -> triplet singlet
	colj.m_j=colj.m_i=0;
	coli.m_i=colij.m_i;
      }
    }
    else {
      if (li->Flav().Strong()) {
        // triplet -> octet triplet
	size_t nc(Flow::Counter());
	coli.m_j=colj.m_i=nc;
	coli.m_i=colij.m_i;
      }
      else {
        // triplet -> singlet triplet
	coli.m_j=coli.m_i=0;
	colj.m_i=colij.m_i;
      }
    }
  }
  else if (lij->Flav().StrongCharge()==-3) {
    if (li->Flav().StrongCharge()==-3) {
      if (lj->Flav().Strong()) {
        // anti-triplet -> anti-triplet octet
	size_t nc(Flow::Counter());
	colj.m_i=coli.m_j=nc;
	colj.m_j=colij.m_j;
      }
      else {
        // anti-triplet -> anti-triplet singlet
	colj.m_j=colj.m_i=0;
	coli.m_j=colij.m_j;
      }
    }
    else {
      if (li->Flav().Strong()) {
        // anti-triplet -> octet anti-triplet
	size_t nc(Flow::Counter());
	coli.m_i=colj.m_j=nc;
	coli.m_j=colij.m_j;
      }
      else {
        // anti-triplet -> singlet anti-triplet
	coli.m_j=coli.m_i=0;
	colj.m_j=colij.m_j;
      }
    }
  }
  else if (lij->Flav().Strong()) {
    if (li->Flav().StrongCharge()==8) {
      if (lj->Flav().StrongCharge()==0) {
        // octet -> octet singlet
	coli.m_i=colij.m_i;
	coli.m_j=colij.m_j;
      }
      else {
        // octet -> octet octet
	size_t nc(Flow::Counter());
	colj.m_i=coli.m_j=nc;
	colj.m_j=colij.m_j;
	coli.m_i=colij.m_i;
      }
    }
    else if (abs(li->Flav().StrongCharge())==3) {
      // octet -> triplet anti-triplet (or vice versa)
      coli.m_i=colij.m_i;
      colj.m_j=colij.m_j;
      if (li->Flav().StrongCharge()<0)
	std::swap<ColorID>(coli,colj);
    }
    else {
      // octet -> singlet octet
      colj.m_i=colij.m_i;
      colj.m_j=colij.m_j;
    }
  }
  else {
    if (abs(li->Flav().StrongCharge())==3) {
      // singlet -> triplet anti-triplet (or vice versa)
      size_t nc(Flow::Counter());
      coli.m_i=colj.m_j=nc;
      if (li->Flav().StrongCharge()<0)
	std::swap<ColorID>(coli,colj);
    }
    else if (li->Flav().StrongCharge()==8) {
      // singlet -> octet octet
      colj.m_i=coli.m_j=Flow::Counter();
      colj.m_j=coli.m_i=Flow::Counter();
    }
    else {
      // singlet -> singlet singlet
      colj.m_i=coli.m_j=0;
      colj.m_j=coli.m_i=0;
    }
  }
  li->SetCol(coli);
  lj->SetCol(colj);
}


size_t Cluster_Amplitude::IdIndex(const size_t &id) const
{
  for (size_t i(0);i<m_legs.size();i++)
    if (m_legs[i]->Id()==id) return i;
  return std::string::npos;
}

namespace ATOOLS {

  class Order_LegID {
  public:
    bool operator()(Cluster_Leg *const a,Cluster_Leg *const b) const
    { return a->Id()<b->Id(); }
  };// end of class Sort_LegID

}

void Cluster_Amplitude::OrderLegs()
{
  std::stable_sort(m_legs.begin(),m_legs.end(),Order_LegID());
}

namespace ATOOLS {

  std::ostream &operator<<
    (std::ostream &ostr,const Cluster_Amplitude &ampl)
  {
    ostr<<ampl.NIn()<<" -> "<<ampl.Legs().size()-ampl.NIn()<<" {\n";
    ostr<<"  \\mu_r = "<<sqrt(ampl.MuR2())
	<<", \\mu_f = "<<sqrt(ampl.MuF2())
	<<", \\mu_q = "<<sqrt(ampl.MuQ2())
	<<", \\mu = "<<sqrt(ampl.Mu2())<<"\n";
    ostr<<"  k_T = "<<sqrt(ampl.KT2())<<", z = "<<ampl.Z()
	<<", phi = "<<ampl.Phi()<<", kin = "<<ampl.Kin()
	<<", K = "<<ampl.LKF()<<"\n";
    ostr<<"  oew = "<<ampl.OrderEW()<<", oqcd = "<<ampl.OrderQCD()
	<<", nlo = "<<ampl.NLO()<<", new = "<<ID(ampl.IdNew())
	<<", ncl = "<<ampl.NewCol()<<", flag = "<<ampl.Flag()<<"\n";
    if (ampl.Decays().size()) {
      std::string ds;
      for (DecayInfo_Vector::const_iterator cit(ampl.Decays().begin());
	   cit!=ampl.Decays().end();++cit) 
        ds+=ToString(**cit)+" ";
      ostr<<"  decs = { "<<ds<<"}\n";
    }
    if (ampl.ColorMap().size()) {
      std::string cs;
      for (CI_Map::const_iterator cit(ampl.ColorMap().begin());
	   cit!=ampl.ColorMap().end();++cit) 
	cs+=ToString(cit->first)+"->"+ToString(cit->second)+" ";
      ostr<<"  cols = { "<<cs<<"}\n";
    }
    for (size_t i(0);i<ampl.Legs().size();++i)
      ostr<<"  "<<*ampl.Legs()[i]<<"\n";
    return ostr<<"}";
  }

  std::ostream &operator<<
    (std::ostream &ostr,const Cluster_Amplitude * ampl)
  {

    ostr<<ampl->NIn()<<" -> "<<ampl->Legs().size()-ampl->NIn()<<" {\n";
    ostr<<"  \\mu_r = "<<sqrt(ampl->MuR2())
        <<", \\mu_f = "<<sqrt(ampl->MuF2())
        <<", \\mu_q = "<<sqrt(ampl->MuQ2())
        <<", \\mu = "<<sqrt(ampl->Mu2())<<"\n";
    ostr<<"  k_T = "<<sqrt(ampl->KT2())<<", z = "<<ampl->Z()
        <<", phi = "<<ampl->Phi()<<", kin = "<<ampl->Kin()
        <<", K = "<<ampl->LKF()<<"\n";
    ostr<<"  oew = "<<ampl->OrderEW()<<", oqcd = "<<ampl->OrderQCD()
        <<", nlo = "<<ampl->NLO()<<", new = "<<ID(ampl->IdNew())
        <<", ncl = "<<ampl->NewCol()<<", flag = "<<ampl->Flag()<<"\n";
    if (ampl->Decays().size()) {
      std::string ds;
      for (DecayInfo_Vector::const_iterator cit(ampl->Decays().begin());
           cit!=ampl->Decays().end();++cit)
        ds+=ToString(**cit)+" ";
      ostr<<"  decs = { "<<ds<<"}\n";
    }
    if (ampl->ColorMap().size()) {
      std::string cs;
      for (CI_Map::const_iterator cit(ampl->ColorMap().begin());
           cit!=ampl->ColorMap().end();++cit)
        cs+=ToString(cit->first)+"->"+ToString(cit->second)+" ";
      ostr<<"  cols = { "<<cs<<"}\n";
    }
    for (size_t i(0);i<ampl->Legs().size();++i)
      ostr<<"  "<<ampl->Legs()[i]<<"\n";
    return ostr<<"}";
  }


}


namespace ATOOLS {
  template class SP(ATOOLS::Cluster_Amplitude);

  template <> Blob_Data<SP(ATOOLS::Cluster_Amplitude)>::~Blob_Data() {}
  template class Blob_Data<SP(ATOOLS::Cluster_Amplitude)>;
  template SP(ATOOLS::Cluster_Amplitude)&
    Blob_Data_Base::Get<SP(ATOOLS::Cluster_Amplitude)>();
}
