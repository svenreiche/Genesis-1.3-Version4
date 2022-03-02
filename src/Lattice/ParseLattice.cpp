#include "ParseLattice.h"

ParseLattice::ParseLattice()
{
}


ParseLattice::~ParseLattice()
{
}


bool ParseLattice::parse(string file, string line, int in_rank) {
    cout << "New Parser started..." << endl;

    rank=in_rank;

    istringstream input;
    string instring;
    ostringstream os;

    ifstream fin(file.c_str(), ios_base::in);

    if (!fin) {
        if (rank == 0) { cout << "*** Error: Cannot open magnetic lattice file: " << file << endl; }
        return false;
    }
    while (getline(fin, instring, '\n')) {
        os << instring << endl;
    }
    fin.close();
    input.str(os.str());


    //------------------------------------------------------
    // step one - coarse parsing of the input deck

    string comstring = "";
    vector<string> content;

    while (getline(input, instring)) {    // read line
        // cout << instring<<endl;
        this->trim(instring);
        if ((!instring.compare(0, 1, "#") || instring.length() < 1)) { continue; } // skip comment and empty rows
        comstring.append(" ");
        comstring.append(instring);  // add all content into one string
    }


    for (int i = 0; i < comstring.size(); i++) { // convert to lower case
        comstring[i] = tolower(comstring[i]);
    }

    size_t pos;
    while ((pos = comstring.find_first_of(";")) != string::npos) {  // split into individual lines
        instring = comstring.substr(0, pos);
        this->trim(instring);
        content.push_back(instring);
        comstring.erase(0, pos + 1);
    }

    this->trim(comstring);
    if ((comstring.length() > 0) && (rank == 0)) {
        cout << "*** Warning: Ignoring incomplete last line in magnetic lattice" << endl;
    }


    // -----------------------------------------------------------------------
    // step two - parse each individual line of the lattice according to the format
    // label: type =(content);
    // e.g.   "QF1: Quadrupole = {L=0.2, k1=0.8};

    // label.clear();
    // type.clear();
    // argument.clear();

    string inlabel, intype, inargument;
    vector<string> args;

    bool error = false;
    rawlat.clear();

    for (int i = 0; i < content.size(); i++) {

        if ((pos = content[i].find_first_of(':')) == string::npos) {
            if (rank == 0) { cout << "*** Error: Invalid Format in lattice file: " << content[i] << endl; }
            error = true;
            continue;
        }
        inlabel = content[i].substr(0, pos);
        this->trim(inlabel);
        content[i].erase(0, pos + 1);
        if ((pos = content[i].find_first_of('=')) == string::npos) {
            if (rank == 0) { cout << "*** Error: Invalid Format in lattice file: " << content[i] << endl; }
            error = true;
            continue;
        }
        intype = content[i].substr(0, pos);
        this->trim(intype);
        content[i].erase(0, pos + 1);

        if ((pos = content[i].find_first_of('{')) == string::npos) {
            if (rank == 0) { cout << "*** Error: Invalid Format in lattice file: " << content[i] << endl; }
            error = true;
            continue;
        }
        content[i].erase(0, pos + 1);
        if ((pos = content[i].find_first_of('}')) == string::npos) {
            if (rank == 0) { cout << "*** Error: Invalid Format in lattice file: " << content[i] << endl; }
            error = true;
            continue;
        }

        content[i].erase(pos, content[i].size());
        inargument = content[i];
        this->trim(inargument);

        args.clear();
        args.push_back("type="+intype.substr(0,4));
        this->parseArguments(args, inargument);
        rawlat[inlabel] = args;
    }
    if (error){ return false; }

    // -------------------------------------------------------------
    // step 3 - resolving all references and convert raw lattice to lattice of type LatElement


    int recursion = 20;
    elements.clear();
    lines.clear();
    for ( auto const &ele : rawlat) {
        string type = ele.second[0].substr(5);
        if (type.compare("line") != 0) {
            error = this->dereference(ele.first, type, recursion);
            if (!error) { return error; }
        } else {
            lines[ele.first]=ele.second;
            reverse(lines[ele.first].begin(),lines[ele.first].end());  // needs flipped order
        }
    }

    for (auto const &ele :elements){
        cout << "Element: " << ele.first << endl;
        for (auto const &par: ele.second){
            cout << "  " << par.first << " = " << par.second << endl;
        }
        cout << endl;
    }


    // -------------------------------------------------------------
    // step 3 - resolving all references and convert raw lattice to lattice of type LatElement

    for (int i=0;i<line.size();i++){
        line[i]=tolower(line[i]);
    }

    if (lines.find(line) == lines.end()){
        if (rank == 0) {cout << "*** Error: Selected beamline " << line << " not define din lattice file." << endl;}
        return false;
    }


    beamline.clear();
    iref = 0;
    recursion = 20;
    return this->unroll(line, recursion);

    // note that beamline is still quite virtual with possible references to variables or sequences
}


bool ParseLattice::unroll(const string line, int recursion){

    // check if maximum level of recursion has been reached
    if (recursion < 0) {
        if (rank == 0) { cout << "*** Error: Too many nested references in beamline" << endl; }
        return false;
    }
    int isave = iref;


    for (auto const &ele : lines[line]) {
        if (ele.find_first_of("type") != string::npos) {
            continue;
        }
        string mele = ele;
        int count = this->checkMultiplier(mele);
        double ref = this->checkReference(mele);
        for (int i = 0; i < count; i++) {
            if (lines.find(mele) == lines.end()) {
                if (elements.find(mele) == elements.end()) {
                    if (rank==0) {cout << "*** Error: Reference to undefined element: " << mele << endl;}
                    return false;
                }
                if (ref < 0) {
                    beamline.push_back({mele, 0, 0, -1});
                } else {
                    beamline.push_back({mele, 0, ref, isave});
                }
                iref++;
            } else {
                bool error = this->unroll(mele, recursion -1);
                if (!error) { return error; }
            }
        }
    }

    return true;
}
double ParseLattice::checkReference(string &element){

    size_t pos=element.find_first_of("@");
    if (pos==string::npos){
        return -1.;
    }
    string num=element.substr(pos+1);
    element.erase(pos);
    this->trim(num);
    this->trim(element);
    return atof(num.c_str());
}

int ParseLattice::checkMultiplier(string &element){

    size_t pos=element.find_first_of("*");
    if (pos==string::npos){
        return 1;
    } else {
        string num=element.substr(0,pos);
        element.erase(0,pos+1);
        this->trim(num);
        this->trim(element);
        return atoi(num.c_str());
    }
}




bool ParseLattice::dereference(const string key, const string type, int recursion) {
     // check if maximum level of recursion has been reached
     if (recursion < 0) {
         if (rank == 0) { cout << "*** Error: Too many nested references" << endl; }
         return false;
     }

     // check if the element has already been defined (by reference)
     if (elements.find(key) != elements.end()) {
         return true;
     }


     // check if element type is known
     if (cele.find(type) == cele.end()) {
         if (rank == 0) {
             cout << "*** Error: Unrecognized element type: " << type << " for element: " << key << endl;
         }
         return false;
     }

     // get default values
     elements[key] = {};
     for (auto const &defval: cele[type]) {
         elements[key][defval.first] = defval.second;
         elements[key]["type"] = type;
     }


     // check if element has reference
     bool hasRef = false;
     string ref;
     for (auto const &arg: rawlat[key]) {
         if (arg.find("ref") != string::npos) {
             size_t pos = arg.find("=");
             if (pos == string::npos) {
                 if (rank == 0) { cout << "*** Error: Invalid argument format for reference: " << arg << endl; }
                 return false;
             }
             hasRef = true;
             ref = arg.substr(pos + 1);
             this->trim(ref);
         }
     }

     // if referred to other element check for type and get data
     if (hasRef) {
         string reftype = rawlat[ref][0].substr(5);
         if (reftype.compare(type) != 0) {
             if (rank == 0) { cout << "*** Error: Mismatched in referred type for element: " << key << endl; }
             return false;
         }
         if (!this->dereference(ref, type, recursion - 1)) { return false; }
         // copy reference into current
         for (auto const &ele: elements[ref]) {
             elements[key][ele.first] = ele.second;
         }
     }

     // parse the argument of the rawlat
     size_t pos;
     for (auto const &arg: rawlat[key]) {
         if ((pos = arg.find_first_of('=')) == string::npos) {
             if (rank == 0) {
                 cout << "*** Error: Invalid format in lattice file for element: " << key << " - " << arg << endl;
             }
             return false;
         }
         string argkey = arg.substr(0, pos);
         string argval = arg.substr(pos + 1);
         this->trim(argkey);
         this->trim(argval);
         if (elements[key].find(argkey) == elements[key].end()) {
             if (rank == 0) {
                 cout << "*** Error: Invalid keyword : " << argkey << " in lattice file for element: " << key << endl;
             }
             return false;
         }
         elements[key][argkey] = argval;
     }
     return true;
 }


void ParseLattice::parseArguments(vector<string> &args, string &inargs)
{
    size_t pos;
    if ((pos = inargs.find_first_of(',')) != string::npos) {
        string residual = inargs.substr(pos + 1);
        this->trim(residual);
        this->parseArguments(args, residual);
        string arg = inargs.substr(0, pos);
        this->trim(arg);
        args.push_back(arg);
    } else {
        this->trim(inargs);
        args.push_back(inargs);
    }
    return;
}











