#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <string.h>
#include <math.h>
#include <ctime>

class AtomType
{
    public:

    std::string type;
    int    atomic_number;
    double sqrt_c6, sqrt_c6_nei;
    double sqrt_c12, sqrt_c12_nei;
    double charge;

    AtomType () {}

    AtomType (std::string type, int atomic_number, double sqrt_c6, double sqrt_c6_nei, double sqrt_c12, double sqrt_c12_nei, double charge)
    {
        this->type = type;
        this->atomic_number = atomic_number;
        this->sqrt_c6 = sqrt_c6;
        this->sqrt_c6_nei = sqrt_c6_nei;
        this->sqrt_c12 = sqrt_c12;
        this->sqrt_c12_nei = sqrt_c12_nei;
        this->charge = charge;
    }

    void Write (FILE* where) 
    {
        char line[2048];
        double c6 = this->sqrt_c6 * this->sqrt_c6;
        double c12 = this->sqrt_c12 * this->sqrt_c12;
        sprintf(line, "%-20s%-10d%-10.3f%-10.3f%-10s%15.7e%15.7e\n", this->type.c_str(), this->atomic_number, 0.0, 0.0, "A", c6, c12);
        fprintf(where, "%s", line);
        return; 
    }
};

class Pair
{
    public:

    AtomType* atom_i;
    AtomType* atom_j;
    double    c6;
    double    c12;

    Pair () {}

    Pair (AtomType* atom_i, AtomType* atom_j)
    {
        this->atom_i = atom_i;
        this->atom_j = atom_j;
        this->c6 = atom_i->sqrt_c6 * atom_j->sqrt_c6;
        this->c12 = atom_i->sqrt_c12 * atom_j->sqrt_c12;
    }

    void Update ()
    {
        this->c6 = atom_i->sqrt_c6 * atom_j->sqrt_c6;
        this->c12 = atom_i->sqrt_c12 * atom_j->sqrt_c12;
    }

    void Write (FILE* where)
    {
        char line[2048];
        sprintf(line, "%-10s%-10s%-10d%-15.7e%-15.7e\n", this->atom_i->type.c_str(), this->atom_j->type.c_str(), 1, this->c6, this->c12);
        fprintf(where, "%s", line);
        return;
    }
};

class PairType
{
    public:

    AtomType* atom_i;
    AtomType* atom_j;
    double    c6;
    double    c12;

    PairType () {}

    PairType (AtomType* atom_i, AtomType* atom_j)
    {
        this->atom_i = atom_i;
        this->atom_j = atom_j;
        this->c6 = atom_i->sqrt_c6_nei * atom_j->sqrt_c6_nei;
        this->c12 = atom_i->sqrt_c12_nei * atom_j->sqrt_c12_nei;
    }

    void Update ()
    {
        this->c6 = atom_i->sqrt_c6_nei * atom_j->sqrt_c6_nei;
        this->c12 = atom_i->sqrt_c12_nei * atom_j->sqrt_c12_nei;
    }

    void Write (FILE* where)
    {
        char line[2048];
        sprintf(line, "%-10s%-10s%-10d%-15.7e%-15.7e\n", this->atom_i->type.c_str(), this->atom_j->type.c_str(), 1, this->c6, this->c12);
        fprintf(where, "%s", line);
        return;
    }
};

class Forcefield
{
    public:

    std::vector<AtomType*> atoms;
    std::vector<Pair*> nonbond_params;
    std::vector<PairType*> pairtypes;

    Forcefield () {}

    void Update ()
    {
        for (int i = 0; i < nonbond_params.size(); i++)
        {
            nonbond_params[i]->Update();
        }
        for (int i = 0; i < pairtypes.size(); i++)
        {
            pairtypes[i]->Update();
        }
        return;
    }

    void UpdatePairs ()
    {
        for (int i = 0; i < nonbond_params.size(); i++)
        {
            nonbond_params[i]->Update();
        }
    }

    void AddAtom (AtomType* atom)
    {
        this->atoms.push_back(atom);
        return;
    }

    AtomType* FindAtom (std::string type)
    {
        for (int i = 0; i < this->atoms.size(); i++)
        {
            if (this->atoms[i]->type == type)
            {
                return this->atoms[i];
            }
        }
        return NULL;
    }

    void ApplyFactorsToAtom (double alpha_c6, double alpha_c12, std::string type)
    {
        // first find the atom
        AtomType* found_atom = this->FindAtom (type);
        if (found_atom != NULL)
        {
            // apply factors
            found_atom->sqrt_c6  *= sqrt(alpha_c6);
            found_atom->sqrt_c12 *= sqrt(alpha_c12);
            // update all pairs of the forcefield, i.e. recalculate c6 and c12 of pairs based on new atomic c6 and c12
            this->Update();
            return;
        }
        // if not found
        else
        {
            fprintf(stderr, "FATAL ERROR: Atom not found.\n");
        }

    }

    // i.e. make pairtypes equal to nonbond_params + self-interactions
    void IgnorePairtypes ()
    {
        for (int i = 0; i < this->atoms.size(); i++)
        {
            this->atoms[i]->sqrt_c6_nei = this->atoms[i]->sqrt_c6;
            this->atoms[i]->sqrt_c12_nei = this->atoms[i]->sqrt_c12;
        }

        // update pairtypes
        for (int i = 0; i < this->pairtypes.size(); i++)
        {
            this->pairtypes[i]->Update();
        }

        return;
    }

    // define    CS_6^{1/2}[i] = F_6[i] * C_6^{1/2}[i]
    // define CS_{12}^{1/2}[i] = F_{12}[i] * C_{12}^{1/2}[i]
    //
    // the vector of factors has to follow the order of definition of the atoms
    // 
    // control option controls which will be altered
    //      = 0     c6
    //      = 1     c12
    //      = 2     c6 + c12
    //
    void SetPairtypesViaFactors (std::vector<double> factors_c6, std::vector<double> factors_c12, int control)
    {
        // check sizes
        if ( (factors_c6.size() != atoms.size()) || (factors_c12.size() != atoms.size()) )
        {
            std::cerr << "error: number of factors for c6 or c12 do not correspond to number of atoms\n";
            exit(0);
        }
        // set factors
        for (int i = 0; i < this->atoms.size(); i++)
        {
            if ( (control==0) || (control==2) )
                this->atoms[i]->sqrt_c6_nei = this->atoms[i]->sqrt_c6 * factors_c6[i];
            if ( (control==1) || (control==2) )
                this->atoms[i]->sqrt_c12_nei = this->atoms[i]->sqrt_c12 * factors_c12[i];
        }
        // update pairtypes
        for (int i = 0; i < this->pairtypes.size(); i++)
        {
            std::cerr << "updating pairtype " << i+1 << "\n";
            std::cerr << "    atom_i = " << pairtypes[i]->atom_i->type << "\n";
            std::cerr << "    atom_j = " << pairtypes[i]->atom_j->type << "\n";
            this->pairtypes[i]->Update();
        }
        return;
    }

    void AddNonbondParam (Pair* pair)
    {
        this->nonbond_params.push_back(pair);
        return;
    }

    void AddPairtype (PairType* pair)
    {
        this->pairtypes.push_back(pair);
        return;
    }

    void FullWrite (FILE* where)
    {
        // header
        // current date/time based on current system
        time_t now = time(0);

        // convert now to string form
        char* dt = ctime(&now);

        fprintf(where, "; This file was generated via ff_lib code on %s; \n", dt);

        fprintf(where,
"[ defaults ]\n"
"; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n"
"  1             1               no              1.0     1.0\n\n");

        // atomtypes
        fprintf(where, "[ atomtypes ]\n");
        for (int i = 0; i < atoms.size(); i++)
        {
            atoms[i]->Write(where);
        }
        fprintf(where, ";\n\n");
        // nonbond_params
        fprintf(where, "[ nonbond_params ]\n");
        for (int i = 0; i < nonbond_params.size(); i++)
        {
            nonbond_params[i]->Write(where);
        }
        fprintf(where, ";\n\n");
        // pairtypes
        fprintf(where, "[ pairtypes ]\n");
        if (this->pairtypes.size() != 0)
        {
            for (int i = 0; i < pairtypes.size(); i++)
            {
                pairtypes[i]->Write(where);
            }
        }
        return;
     }

    void Write (FILE* where)
    {
        // header
        // current date/time based on current system
        time_t now = time(0);

        // convert now to string form
        char* dt = ctime(&now);

        fprintf(where, "; This file was generated via ff_lib code on %s; \n", dt);

        // atomtypes
        fprintf(where, "[ atomtypes ]\n");
        for (int i = 0; i < atoms.size(); i++)
        {
            atoms[i]->Write(where);
        }
        fprintf(where, ";\n\n");
        // nonbond_params
        fprintf(where, "[ nonbond_params ]\n");
        for (int i = 0; i < nonbond_params.size(); i++)
        {
            nonbond_params[i]->Write(where);
        }
        fprintf(where, ";\n\n");
        // pairtypes
        fprintf(where, "[ pairtypes ]\n");
        if (this->pairtypes.size() != 0)
        {
            for (int i = 0; i < pairtypes.size(); i++)
            {
                pairtypes[i]->Write(where);
            }
        }
        return;
     }

    void Write (std::string fn)
    {
        FILE* fp = fopen(fn.c_str(), "w");
        this->Write(fp);
        fclose(fp);
    }

    // FullWrite includes the force-field headers, but not the bonded parameters.
    // It is assumed that the bonded parameters are specified in the topology itself.
    void FullWrite (std::string fn)
    {
        FILE* fp = fopen(fn.c_str(), "w");
        this->FullWrite(fp);
        fclose(fp);
    }

    void ReadAtomsFromInput (std::string fn)
    {
        char buffer[2048];
        char type[4];
        int  an;
        double sqrt_c6, sqrt_c6_nei, sqrt_c12, sqrt_c12_nei;
        AtomType* new_atom;
        FILE* fp = fopen(fn.c_str(), "r");
        while (fgets(buffer, 2048, fp))
        {
            if (buffer[0] == '#')
                continue;
            else
            {
                sscanf(buffer, "%s%d%lf%lf%lf%lf", type, &an, &sqrt_c6, &sqrt_c12, &sqrt_c6_nei, &sqrt_c12_nei); 
                new_atom = new AtomType (type, an, sqrt_c6, sqrt_c6_nei, sqrt_c12, sqrt_c12_nei, 0.00);
                this->AddAtom(new_atom);
            }
        }
        fclose(fp);
        return;
    }

    void GeneratePairsFromAtoms ()
    {
        Pair* new_pair;
        PairType* new_pairt;
        for (int i = 0; i < this->atoms.size(); i++)
        {
           for (int j = i; j < this->atoms.size(); j++)
           {
               if (j > i)
               {
                   // nonbond -- only needed for j != i, otherwise it is defined in atomtypes
                   new_pair = new Pair(this->atoms[i], this->atoms[j]);
                   this->nonbond_params.push_back(new_pair);
               }
               // pairtypes
               new_pairt = new PairType(this->atoms[i], this->atoms[j]);
               this->pairtypes.push_back(new_pairt);
           }
        }
    }

};
