package ProstateCancer;

import HAL.GridsAndAgents.AgentGrid2D;
import HAL.GridsAndAgents.AgentSQ2D;
import HAL.GridsAndAgents.PDEGrid2D;
import HAL.Gui.GridWindow;
import HAL.Rand;
import HAL.Tools.FileIO;
import HAL.Util;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import static HAL.Util.*;
import static java.lang.Integer.parseInt;

class ProstateCancerCell extends AgentSQ2D<ProstateCancerGrid> { //the agents
    public int color;
    public void Init(double f0){//Method to assign random color to each cell
        if (G.rng.Double()<f0) {
            this.color = RED;
        }else{
            this.color = GREEN;
        }
    }
    public void StepCel(){
        if (this.color ==GREEN){
            if (G.rng.Double() < (G.divProb + G.dieProb)/(1-G.moveProb)) {//Check if any event occur (G.divProb*divProbS + G.dieProb)
                if (G.rng.Double() < G.divProb / (G.divProb + G.dieProb)) { //G.divProb*divProbS / (G.divProb*divProbS + G.dieProb)
                    int options = G.MapEmptyHoodK(G.divHood,Xsq(),Ysq());//Checking the available potions in the neighbourhood
                    if (options > 0) {
                        if (G.rng.Double() > G.dieProbD * G.GF.Get(Isq())) {//Decide whether to dispose a cell
                            //Dispose();  G.NewAgentSQ(Isq()).color = this.color;
                            int loc = G.rng.Int(options);
                            if(G.PopAt(G.divHood[loc])<G.K ){
                                G.NewAgentSQ(G.divHood[loc]).color=this.color;
                            }
                        } else {
                            Dispose();
                        }
                    }
                } else {
                    Dispose();
                }
            }
        }else {
            if (G.rng.Double() < (G.divProbR + G.dieProb)/(1-G.moveProb)) {//Check if any event occur//(G.divProb*divProbR + G.dieProb)
                if (G.rng.Double() < G.divProbR / (G.divProbR + G.dieProb)) { //G.divProb*divProbR / (G.divProb*divProbR + G.dieProb)
                    int options = G.MapEmptyHoodK(G.divHood,Xsq(),Ysq());//Checking the available potions in the neighbourhood
                    if (options > 0) {
                        if (G.rng.Double() >= 0.0 * G.GF.Get(Isq())) {//Decide whether to dispose a cell
                            //Dispose(); G.NewAgentSQ(Isq()).color = this.color;
                            int loc = G.rng.Int(options);
                            if(G.PopAt(G.divHood[loc])<G.K ){
                                G.NewAgentSQ(G.divHood[loc]).color=this.color;
                            }
                        } else {
                            Dispose();
                        }
                    }
                } else {
                    Dispose();
                }
            }
        }
    }
    public void StepCelFibro(){
          if (this.color ==GREEN){
            if (G.rng.Double() < (G.divProbFibro + G.dieProb)/(1-G.moveProb)) {//Check if any event occur (G.divProb*divProbS + G.dieProb)
                if (G.rng.Double() < G.divProbFibro / (G.divProbFibro + G.dieProb)) { //G.divProb*divProbS / (G.divProb*divProbS + G.dieProb)
                    int options = G.MapEmptyHoodK(G.divHood,Xsq(),Ysq());//Checking the available potions in the neighbourhood
                    if (options > 0) {
                        if (G.rng.Double() > G.dieProbD * G.GF.Get(Isq())) {//Decide whether to dispose a cell
                            //Dispose();  G.NewAgentSQ(Isq()).color = this.color;
                            int loc = G.rng.Int(options);
                            if(G.PopAt(G.divHood[loc])<G.K ){
                                G.NewAgentSQ(G.divHood[loc]).color=this.color;
                            }
                        } else {
                            Dispose();
                        }
                    }
                } else {
                    Dispose();
                }
            }
        }else {
            if (G.rng.Double() < (G.divProbRFibro + G.dieProb)/(1-G.moveProb)) {//Check if any event occur//(G.divProb*divProbR + G.dieProb)
                if (G.rng.Double() < G.divProbRFibro / (G.divProbRFibro + G.dieProb)) { //G.divProb*divProbR / (G.divProb*divProbR + G.dieProb)
                    int options = G.MapEmptyHoodK(G.divHood,Xsq(),Ysq());//Checking the available potions in the neighbourhood
                    if (options > 0) {
                        if (G.rng.Double() >= 0.0 * G.GF.Get(Isq())) {//Decide whether to dispose a cell
                            //Dispose(); G.NewAgentSQ(Isq()).color = this.color;
                            int loc = G.rng.Int(options);
                            if(G.PopAt(G.divHood[loc])<G.K ){
                                G.NewAgentSQ(G.divHood[loc]).color=this.color;
                            }
                        } else {
                            Dispose();
                        }
                    }
                } else {
                    Dispose();
                }
            }
        }
    }

    public void moveCel(){
        if (this.IsAlive()){
            //displaces the particle, particle will not move if movement would cause it to move out of the grid
            int options = G.MapEmptyHoodK(G.divHood,Xsq(),Ysq());
            if (options>0){
                MoveSQ(G.divHood[G.rng.Int(options)]);
            }
        }
    }
}

public class ProstateCancerGrid extends AgentGrid2D<ProstateCancerCell> {//Grid to house the agents
    double dT = 1;
    boolean dec; boolean AT = false; int N0, rCount, rs, K; double FibroImpact = 200;
    double divProb = 0.027;//base growth rate
    double divProbR = (1-0.3)*divProb;
    double dieProb = 0.3*divProb;//constant death rate
    double dieProbD = 0.75;
    double moveProb = 1*divProb;
    double divProbFibro = FibroImpact*divProb/100;
    double divProbRFibro = FibroImpact*divProbR/100;
    public PDEGrid2D GF;//The PDEGrid2D to model the growth factor(GF) diffusion
    double initPA; double PA;  double growthPAS = 1; double growthPAR = 1; double decayPA = 0.5;
    Rand rng = new Rand();// Define rng as a random number generator object
    int[]divHood= Util.VonNeumannHood(true);//RectangleHood(false,2,2);//
    // //counting the neighborhood
    public ProstateCancerGrid(int x, int y) {//Construct grid to house the agents
        super(x, y, ProstateCancerCell.class);
        GF = new PDEGrid2D(x,y);//Defining PDE grid to solve the PDE model
    }
    public int MapEmptyHoodK(int[] hood, int centerX, int centerY) {
        return MapHood(hood, centerX, centerY, (i, x, y) -> PopAt(i) < K );
    }
    public void StepCells(){
        GF.Update();//Update the PDE solution
        ShuffleAgents(rng);//Shuffle the order of simulation
        for(ProstateCancerCell cell: this) {//Run the ABM for each cell
            if (rng.Double()>this.moveProb) {
                if //((cell.Xsq()>12&&cell.Ysq()>31&&cell.Xsq()<23&&cell.Ysq()<42)||(cell.Xsq()>54&&cell.Ysq()>47&&cell.Xsq()<65&&cell.Ysq()<68)||(cell.Xsq()>88&&cell.Ysq()>79&&cell.Xsq()<99&&cell.Ysq()<90)||(cell.Xsq()>0&&cell.Ysq()>72&&cell.Xsq()<11&&cell.Ysq()<83)||(cell.Xsq()>7&&cell.Ysq()>81&&cell.Xsq()<18&&cell.Ysq()<92)||(cell.Xsq()>38&&cell.Ysq()>63&&cell.Xsq()<49&&cell.Ysq()<74)||(cell.Xsq()>7&&cell.Ysq()>41&&cell.Xsq()<18&&cell.Ysq()<32)||(cell.Xsq()>47&&cell.Ysq()>34&&cell.Xsq()<58&&cell.Ysq()<44)||(cell.Xsq()>63&&cell.Ysq()>19&&cell.Xsq()<74&&cell.Ysq()<30)||(cell.Xsq()>45&&cell.Ysq()>19&&cell.Xsq()<56&&cell.Ysq()<30)) {
                  // (cell.Xsq()>15&&cell.Ysq()>34&&cell.Xsq()<49&&cell.Ysq()<67) {//SQO
                 //  (cell.Xsq()>0&&cell.Ysq()>34&&cell.Xsq()<33&&cell.Ysq()<67) {//SQS
                // ((cell.Xsq()>20&&cell.Ysq()>24&&cell.Xsq()<31&&cell.Ysq()<81)||(cell.Xsq()>30&&cell.Ysq()>70&&cell.Xsq()<77&&cell.Ysq()<81)) {//OU1
                  //  ((cell.Xsq()>10&&cell.Ysq()>34&&cell.Xsq()<21&&cell.Ysq()<91)||(cell.Xsq()>20&&cell.Ysq()>80&&cell.Xsq()<67&&cell.Ysq()<91)) {//OU2
                   // ((cell.Xsq()>22&&cell.Ysq()>22&&cell.Xsq()<73&&cell.Ysq()<28)||(cell.Xsq()>72&&cell.Ysq()>22&&cell.Xsq()<78&&cell.Ysq()<73)||(cell.Xsq()>27&&cell.Ysq()>72&&cell.Xsq()<78&&cell.Ysq()<78)||(cell.Xsq()>22&&cell.Ysq()>27&&cell.Xsq()<28&&cell.Ysq()<78)) { //FSQ2
                    //((cell.Xsq()>32&&cell.Ysq()>32&&cell.Xsq()<58&&cell.Ysq()<43)||(cell.Xsq()>57&&cell.Ysq()>32&&cell.Xsq()<68&&cell.Ysq()<58)||(cell.Xsq()>42&&cell.Ysq()>57&&cell.Xsq()<68&&cell.Ysq()<68)||(cell.Xsq()>32&&cell.Ysq()>42&&cell.Xsq()<43&&cell.Ysq()<68)) {//
                    (cell.Xsq()>34&&cell.Ysq()>34&&cell.Xsq()<67&&cell.Ysq()<67) { //FC
                    cell.StepCelFibro();
                }else{
                    cell.StepCel();
                }
            }else{
                cell.moveCel();
            }
        }
    }

    public void DrawModel(GridWindow win, int iPatch){//Drawing model
        for (int x = 0; x < xDim; x++) {//loop through all x
            for (int y = 0; y < yDim; y++) {//loop through all y
                if(GetAgent(x,y)!=null){
                    win.SetPix(x + iPatch * xDim,y,GetAgent(x,y).color);//win.SetPix(xDim+x,y,color);
                }else{
                    win.SetPix(x + iPatch * xDim,y, BLACK);
                }
            }

        }
    }
    public static void main(String[]args) {
        int x = 100, y = 100, timesteps = 2000, nIter = 30;
        double f0 = 0.1; int patchCount = 1; int redLevel = 2;

        GridWindow win = new GridWindow("Sensitive Cell(Green), Resistant cell(Red)", x*patchCount, y, 3);//Define the display window
        for (int iIter = 0; iIter < nIter; iIter++) {
            ProstateCancerGrid[] model = new ProstateCancerGrid[patchCount];//Define the model

            for (int k = 1; k < 2; k++) {
                for (int ip = 0; ip < patchCount; ip++) {
                    model[ip] = new ProstateCancerGrid(x, y);
                    model[ip].N0 = 5000;

                    String path = "C:\\HAL-master\\ProstateCancer\\clumpInput_010.csv";
                    int lineCount = 0;
                    model[ip].rCount = 0;
                    String line = "";
                    try {
                        BufferedReader br = new BufferedReader(new FileReader(path));
                        while ((line = br.readLine()) != null) {
                            String[] values = line.split(",");
                            for (int i = 0; i < values.length; i++) {
                                System.out.println(values[i]);
                                if (parseInt(values[i]) == 1) {
                                    model[ip].NewAgentSQ(lineCount * x + i).color = GREEN;
                                    System.out.println("G");
                                }
                                if (parseInt(values[i]) == 2) {
                                    model[ip].NewAgentSQ(lineCount * x + i).color = RED;
                                    model[ip].rCount++;
                                    System.out.println("R");
                                }
                            }
                            lineCount++;
                        }
                    } catch (FileNotFoundException e) {
                        e.printStackTrace();
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
                model[0].initPA = (int) 2 * model[0].N0;
                model[0].dec = true;
                model[0].K = k;
                String filename = "CT_K_" + k + "_m" + model[0].moveProb + "_rho" + model[0].growthPAS + redLevel + model[0].divProbR + "_FC_" + model[0].FibroImpact+ "clumpInput_010";

                FileIO popsOut = new FileIO(filename + "_" + iIter + ".csv", "w");
                popsOut.Write("divProbS" + "," + "divProbR" + "," + "dieProb" + "," + "moveProb" + "," + "initPA" + "," + "growthPAS" + "," + "growthPAR" + "," + "decayPA" + "," + "f0" + "," + "redLevel" + "\n");
                popsOut.Write(model[0].divProb + "," + model[0].divProbR + "," + model[0].dieProb + "," + model[0].moveProb + "," + model[0].initPA + "," + model[0].growthPAS + "," + model[0].growthPAR + "," + model[0].decayPA + "," + f0 + "," + redLevel + "\n");
                popsOut.Write(model[0].N0 + "," + model[0].rCount + "," + model[0].initPA + "\n");
                popsOut.Write("Dose(1)" + "," + "Sensitive(1)" + "," + "Resistant(1)" + "," + "Total(1)" + "," + "Biomarker(1)" + "\n");

                double tPA = model[0].initPA;

                for (int i = 0; i < timesteps; i++) {

                    if (model[0].AT == false) {//Continuous Therapy
                        System.out.println(k + "" + model[0].dec + "0:" + i);
                        model[0].GF.SetAll(1);
                    } else {//Adaptive Therapy
                        for (int ip = 0; ip < patchCount; ip++) {
                            if (model[ip].Pop() > model[ip].N0 / redLevel && model[ip].dec) {
                                model[ip].GF.SetAll(1);
                            } else if (model[ip].Pop() <= model[ip].N0 / redLevel) {
                                model[ip].dec = false;
                                model[ip].GF.SetAll(0);
                            } else if (model[ip].Pop() > model[ip].N0) {
                                model[ip].dec = true;
                            } else {
                                model[ip].GF.SetAll(0);
                            }
                            System.out.println(k + "" + model[ip].dec + ip + ":" + i);
                        }
                    }
                    tPA = 0.0;
                    for (int ip = 0; ip < patchCount; ip++) {//Simulating the model over the time span
                        win.TickPause(100);
                        //model step
                        model[ip].StepCells();//Execute the model for each time step
                        //draw
                        model[ip].DrawModel(win, ip);//Show the output
                        if ((i + 1) % 1 == 0) {
                            FileIO distOut = new FileIO(filename + "_time_" + i + "_" + iIter + ".csv", "w");
                            for (int iy = 0; iy < model[ip].yDim; iy++) {
                                for (int ix = 0; ix < model[ip].xDim; ix++) {
                                    if (model[ip].GetAgent(ix, iy) != null) {
                                        if(model[ip].PopAt(ix,iy)==1){
                                            if (model[ip].GetAgent(ix, iy).color == GREEN) {
                                                distOut.Write("1");
                                            } else {
                                                distOut.Write("2");
                                            }
                                        } else{
                                            model[ip].rs = 0;
                                            for (ProstateCancerCell cell : model[ip].IterAgents(ix, iy)) {//iterates over each agent in the grid point to count the number of all red cells in the grid
                                                if (cell.color == RED) {
                                                    model[ip].rs++;
                                                }
                                            }
                                            if(model[ip].rs==0){
                                                distOut.Write("3");
                                            }else if(model[ip].rs==1){
                                                distOut.Write("4");
                                            }else{
                                                distOut.Write("5");
                                            }
                                        }
                                    }else {
                                        distOut.Write("0");
                                    }
                                    if (ix < model[ip].xDim - 1) distOut.Write(",");
                                }
                                distOut.Write("\n");//System.out.println(iy);
                            }
                            distOut.Close();
                        }
                        model[ip].rCount = 0;//Counting cell population
                        for (int ix = 0; ix < model[ip].xDim; ix++) {
                            for (int iy = 0; iy < model[ip].yDim; iy++) {
                                if (model[ip].GetAgent(ix, iy) != null) {
                                    for (ProstateCancerCell cell : model[ip].IterAgents(ix, iy)) {//iterates over each agent in the grid point to count the number of all red cells in the grid
                                        if (cell.color == RED) {
                                            model[ip].rCount++;
                                        }
                                    }
                                }
                            }
                        }//Counting cell population ends
                        //Counting Distance
                        model[ip].PA = model[ip].PA + model[ip].growthPAS * (model[ip].Pop() - model[ip].rCount) + model[ip].growthPAR * model[ip].rCount - model[ip].decayPA * model[ip].PA;
                        tPA = tPA + model[ip].PA;
                        popsOut.Write(model[ip].GF.GetAvg() + "," + (model[ip].Pop() - model[ip].rCount) + "," + model[ip].rCount + "," + model[ip].Pop() + "," + model[ip].PA);
                    }
                    popsOut.Write("\n");
                }
                popsOut.Close();
            }
        }
        win.Close();
    }

}
