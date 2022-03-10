#include "ParseLattice.h"


ParseLattice::ParseLattice()
{
}


ParseLattice::~ParseLattice()
{
}

void ParseLattice::generateLattice(double dz, map<string,vector<double>> &lattice) {


    // clear first old lattice
    for (auto &ele : lattice){
        ele.second.clear();
    }

    double z=0;
    double l;
    vector<double> zloc;
    vector<pair<double,double>> zfine;
    zloc.push_back(0);

    // generate first the unordered lattice of the starting and end point of all elements
    for (auto& ele: beamline){
        // check if beam line has absolute position
        if (ele.ref >=0) {
            z = beamline[ele.ref].z + ele.zoff;
            zloc.push_back(z);
        }

        l = 0;
        if (elements[ele.key]["type"] == "undu") {
            l = atof(elements[ele.key]["lambdau"].c_str()) * atoi(elements[ele.key]["nwig"].c_str());
            zfine.push_back({z,l});
        }
        else {
            if (elements[ele.key].find("l") != elements[ele.key].end()){
                l = atof(elements[ele.key]["l"].c_str());
            }
        }
        ele.z = z;
        ele.l = l;
        z+=l;
        zloc.push_back(z);
    }

    // sort the entries
    sort(zloc.begin(),zloc.end());

    // fill the final lattice grid and eliminate zero elements (e.g. markers)
    double zlast=-1.;
    for (auto const &ele : zloc){
        if (ele > zlast) {
            lattice["z"].push_back(ele);
            zlast=ele;
        }
    }

    // search for undulators and add additional lattice points
    for (auto const & ele :zfine) {
        cout << "Undulator from " << ele.first << " to " << ele.first + ele.second << endl;
        auto first = find(zloc.begin(), zloc.end(), ele.first);
        auto last = find(zloc.begin(), zloc.end(), ele.first + ele.second);
        for (auto iter = first; iter < last; iter++) {
            auto delz=*(iter+1)-*(iter);
            if (delz > 0) {
                int nfine = static_cast<int>(round(delz / dz));
                if (nfine < 1) { nfine = 1; }
                auto dzfine = delz / static_cast<double>(nfine);
                cout << " Distances: " << delz << " with step size (mm) " << dzfine * 1e3 << " for " << nfine <<
                    " steps (Reference: " << dz * 1e3 << " )" << endl;
                for (int i =1; i < nfine ; i ++){
                    lattice["z"].push_back(*iter+dzfine*static_cast<double>(i) );
                }
            }
        }
    }
    sort(lattice["z"].begin(),lattice["z"].end());

    // allocate space for other lattice objects and initialize them to zero
    int nlat = lattice["z"].size();
    cout << "Lattice size is " << nlat << endl;

    for (auto &ele: lattice){
        if (ele.first != "z") {
            ele.second.resize(nlat-1);
            for (auto &val: ele.second) {
                val = 0;
            }
        }
    }
    cout << " allocate: " << lattice["aw"].size() << endl;
    map<string,double> args;

    // fill lattice with all other elements
    for (auto const &ele: beamline){
        double z0 = ele.z;
        double z1 = z0 + ele.l;
        double temp = 0;
        string type = this->convertArguments(ele.key, args);
        cout << "Parsing " << ele.key << " of type " << type << endl;
        for (auto const & arg : args) {
            cout << "   " << arg.first << " : " << arg.second << endl;
        }
        if (type=="chic"){
            temp = this->calculateChicAngle(args["delay"],args["lb"],args["ld"]);
        }
        int first = static_cast<int> (find(lattice["z"].begin(), lattice["z"].end(), z0)- lattice["z"].begin());
        int last = static_cast<int> (find(lattice["z"].begin(), lattice["z"].end(), z1)- lattice["z"].begin());
        cout << " Index range: " << first << " to "  << last <<  " (" << z0 << " to " << z1 << ")" << endl;
        for (int idx = first; idx < last; idx++){
            if (type=="undu") {
                lattice["aw"][idx]=args["aw"];
                lattice["kx"][idx]=args["kx"];
                lattice["ky"][idx]=args["ky"];
                lattice["ax"][idx]=args["ax"];
                lattice["ay"][idx]=args["ay"];
                lattice["gradx"][idx]=args["gradx"];
                lattice["grady"][idx]=args["grady"];
                double twopi = 4*asin(1.);
                if (args["lambdau"] > 0) { lattice["ku"][idx]=twopi/args["lambdau"];}
                lattice["helical"][idx]=args["helical"];
            } else if (type=="quad") {
                lattice["qf"][idx]=args["k1"];
                lattice["qx"][idx]=args["dx"];
                lattice["qy"][idx]=args["dy"];
            } else if (type=="corr") {
                lattice["cx"][idx]=args["cx"];
                lattice["cy"][idx]=args["cy"];
            } else if (type=="phas") {
                lattice["phaseshift"][idx]=args["phi"];
            } else if (type=="chic"){
                lattice["chic_angle"][idx] = temp;
                lattice["chic_lb"][idx]=args["lb"];
                lattice["chic_ld"][idx]=args["ld"];
                lattice["chic_lt"][idx]=args["l"];
            } else if (type=="mark") {
                lattice["marker"][idx]=(args["dumpfield"] > 0 ? 1 : 0)
                                        +(args["dumpbeam"]> 0 ? 2 : 0 )
                                        +(args["sort"]>0 ? 4 : 0)
                                        +(args["stop"]>0 ? 8 : 0);
            }
        }
    }

    for (int idx =0; idx < nlat-1; idx++){
        lattice["dz"][idx]=lattice["z"][idx+1]-lattice["z"][idx];
    }


    return;
}

double ParseLattice::calculateChicAngle(double delay0, double lb, double ld){
    if (delay0==0){
        return 0.;
    }
    double delay=fabs(delay);
    double tmin=0;
    double tmax=asin(1)-0.001;
    bool converged=false;
    double theta,d;
    while (!converged) {
        theta = 0.5 * (tmax + tmin);
        d = 4 * lb * (theta / sin(theta) - 1) + 2 * ld * (1 / cos(theta) - 1);
        if (d > delay) {
            tmax = theta;
        } else {
            tmin = theta;
        }
        if (fabs(delay - d) < 1e-15) { converged = true; }
    }
    return delay;
}



string ParseLattice::convertArguments(string key, map<string, double> &args){
    args.clear();
    // here comes the check for variables or sequences
    for (auto const &arg : elements[key]){
        if (arg.first != "type") { args[arg.first]=atof(arg.second.c_str());}
    }
    return elements[key]["type"];
}


bool ParseLattice::parse(string file, string line, int in_rank) {


    rank=in_rank;

    if (rank==0) {  cout << "New Parser started..." << endl; }
    
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
    if (rank == 0) {
    for (auto const &ele :elements){
        cout << "Element: " << ele.first << endl;
        for (auto const &par: ele.second){
            cout << "  " << par.first << " = " << par.second << endl;
        }
        cout << endl;
    }
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
                    beamline.push_back({mele, 0, 0,0, -1});                // z and l are define latter
                } else {
                    beamline.push_back({mele, 0, 0, ref, isave});
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











