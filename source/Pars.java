class Pars {
	/********************************************************************************************
	* MAIN OPTIONS    **************************************************************************
	********************************************************************************************/

	//timing
	static int frameTime = 1; // frame time in minutes
    static int movTime = 24*60;//rate of movie output - frames
    static int[] collectTimesB = new int[]{5*24*60,10*24*60,17*24*60};//times to collect bulk data
    static int[] collectTimesI = new int[]{2*24*60,10*24*60,17*24*60};//times to collect individual cell data
    static int[] dataTime = new int[]{13*60,25*60,24*60};//duration of data collection
    static int txSt = 14*60*24;//treatment start
    static int txDur = 28*60*24;//28*60*24;//treatment duration
	static int[] collectTimesBTx = new int[]{txSt-5*24*60,txSt-24*60,txSt+txDur-24*60};//times to collect during Txs

	//tracks
	static boolean collect = false;//start metrics collection at false
	static int maxTracks = 200;//number of cells to track
	static float rimPercent = 0.25f;//threshold density of carrying capacity to track cells

    //field stuff
    static float diffGW = 0.66f;//defines the fraction for carrying capacity in white matter
    static int maxMesh = 40; //mesh carrying capacity - does not adjust with size

	//layout
    static final int MAX_CELLS = 1000000; // Max number allowed cells, will quit ow
	static final int sizeW = 833;//width of domain in pixels
	static final int sizeH = 573;//height of domain in pixels
	static float sizeMicm = 10000;//height of domain in microns
	static int centerX = Pars.sizeW/2+143;//where to start tumor - x
	static int centerY = Pars.sizeH/2-129;//where to start tumor - y
	static int disp = (int) (3.f*(Pars.sizeMicm/5000.f));//initial dispersion of cells
	static int cellScale = 6;//scale factor for displaying cells
	static int concScale = 5;//scale factor for concentration intensity display

	static float micMToPix = sizeH/(sizeMicm); //microns to pixels conversion
    private static float diameter = 25;//cell diameter in microns
	static float rad = diameter*micMToPix/2; //radius of cell; microns to pixels
	static float hexSide = rad; //percent of cell radius for hexagon side
	static float hexDiag = (float) ((Math.sqrt(3))*hexSide/2);//the line from middle to middle of side
	static int numPointsW = (int) Math.floor((sizeW)/(3*hexSide));//number of points in hex grid - x
	static int numPointsH = (int) Math.floor((sizeH)/(hexDiag));//number of points in hex grid - y
	static double concCutoff = .0005;//concentration of gF cutoff to consider for activation

	//white/gray field information
	static int[][] fArray = new int[sizeW][sizeH];//white matter array
	static int[][] f2Array = new int[sizeW][sizeH];//no matter array

	//region of interest for histology and Ki67 staining
	static int histSizeX = (int) (500.f*Pars.micMToPix);//width
	static int histSizeY = (int) (400.f*Pars.micMToPix);//height
	static int histScale = 10;//scaling of region for graphics
	static int histX=(int) (2.f*Pars.sizeW/3-10);//x location
	static int histY=(int) (1.f*Pars.sizeH/3-80);//y location

	//persistence arrays
    static int sizeSArray = 339;
	static int sizeGArray = 450;
    static int[] stopArray = new int[sizeSArray];//stop times distribution
    static int[] goArray = new int[sizeGArray];//go times distribution

	//conversions
	static float dCConv = frameTime*(10000.f*micMToPix*10000.f*micMToPix)/(24.f*60.f);//convert cm^2/day to um^2/min
	static int divConv = (int) (60.f/frameTime);//converts division times from hours to frames
	static float speedConv = Pars.micMToPix*Pars.frameTime/60.f;//converts speeds from microns/h to pixels/frame
	static float consRConv = frameTime/(24f*60);//convert from per day to per frame

	/*******************
	* VALUES SET FROM INPUT
	******************/
	static String outFile; //output file name
	static int tumN;
	static boolean movie;//record movie?
	static boolean tracks;//record tracks?
	static boolean phenos;//record phenotypes?
	static int txType;//treatment type

	//parameters changed with PDGF
	static float Dc;//diffusion coefficient
	static float startCons;//initial PDGF bolus
	static float pDecay;//PDGF decay rate
	static float secR;//PDGF secretion rate
	static float consR;//PDGF consumption rate
	static float halfMaxP;//concentration of PDGF at half max activity for proliferation
	static float halfMaxM;//concentration of PDGF at half max activity for migration
    static float boost;//concentration of autocrine boost for infected cells
    static float recP;//factor to adjust half max for recruited cell proliferation
    static float recM;//factor to adjust half max for recruited cell migration

    //angle distribution
    static float theEr;//std dev angle white matter

	//SINGLE CELL PARAMETERS

    //proliferation distributions - times in hours
	static int dMax = 200;//max div hours
	private static int dMin = 20;//min div hours
	static int dBins = 20;//number of bins for div
	static int div;//division rate
	static int divEr;//division std dev

	//convert from hours to frames
	static int divMax = dMax*divConv;
	static int divMin = dMin*divConv;
	static int divMean = div*divConv;
	static int divError = divEr*divConv;

	//speed - speeds in microns per hour
	private static float sMax = 180f;//max speed microns/h
	private static float sMin = 0.000f;//min speed microns/h
	static int sBins = 20;//number bins speed
	static float sp;//migration speed
	static float spEr;//migration std dev

	//convert from per hour to per frame
	static float spMax = sMax*speedConv;
	static float spMin = sMin*speedConv;
	static float spMean = sp*speedConv;
	static float spError = spEr*speedConv;
	 	    
}
