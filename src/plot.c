/*****************************************************************************************************************************/
/**															    **/
/**															    **/
/**	Enough. It can't be that you need twenty libraries and several Mbytes to make a bloody xy plot ...		    **/
/**															    **/
/**															    **/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
/**															    **/
/**															    **/
/**	    If you think that this program is unreadable, it is only because it is indeed unreadable.                       **/
/**															    **/
/**															    **/
/*****************************************************************************************************************************/

/*****************************************************************************************************************************/
/**															    **/
/**	Nov 27th, 2012:  - Started writing down changes ... :-)     							    **/
/**	                 - Flag '-hc' prepares commulative histograms    						    **/
/**	Dec 21st, 2012:  - Kolmogorov-Zurbenko filtering added for xy plots (key S)					    **/
/**	Dec 19th, 2015:  - When three columns and the last integer, use it as color coding for dot plots		    **/
/** 	Nov 21st, 2016:  - Freedman-Diaconis rule for estimating number of histogram bins				    **/
/**	Jan 31st, 2017:  - Tiny bug in histograms, plus proper histogram plotting.					    **/
/**	Mar  9th, 2017:  - Flag '-xydy' to allow plotting of standard deviations					    **/
/**	Mar 22nd, 2017:  - When doing a scatter plot, automatically produce density and log-density distribution matrices   **/
/**	Mar 23rd, 2017:  - Flag '-f' to produce density maps on a fine grid (as obtained from Freedman-Diaconis).	    **/
/**	Apr  2nd, 2017:  - Re-nice labeling. Much better.				 				    **/
/**	Apr 20th, 2017:  - Small bug with gl2ppm when drawing histograms and standard deviations			    **/
/**	Jun  9th, 2018:  - Flag "-A" to autoscale second plot to first when overlaying two plots			    **/
/**	Jun 12th, 2018:  - KZ filter now uses Freedman-Diaconis rule (and made much faster)				    **/
/**	Jul 17th, 2018:  - Key 'T' performs and plots the DFT of data (naive, not FFT)      			            **/
/**			 - Scroll with arrow keys, zoom-in/out with page-up/down keys					    **/
/**	Oct  3rd, 2018:  - Changed default for labels, fixed a couple of small bugs					    **/
/**			 - Calculate and print linear correlation coefficient for scatter plots				    **/
/**	Feb  7th, 2019:  - Fix (but only partly) an annoying bug at high zoom-in levels. Not a proper/complete fix. 	    **/
/**			 - Add '-d' flag to allow selection of concatenated data sets (for Qs mainly) 		            **/
/**															    **/
/**															    **/
/**				    SWITCHED TO GITHUB, THE REST FROM THE COMMITS ...                           	    **/
/**															    **/
/**															    **/
/*****************************************************************************************************************************/

#include <Ygl.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <limits.h>
#include <ctype.h>                     
#include <time.h>
#include <omp.h>

#define VERSION "1.201"

#define	NO	0
#define	YES	1
#define	SMALLM	1000.0
#define FREE_ARG char*
#define NR_END  1

#define         MAXRAND         ((float)(INT_MAX)+1.0)
#define         getrand         ((float)(random())/MAXRAND)

extern void  arcxf              (  Coord,  Coord,  Coord,  Coord, Angle, Angle );
void 	sort( int n, float *ra);
void    welford(int N, float *data);
void    myexit( int code );



/* For color calculation */
float	  RED_1 = 0.0;
float	GREEN_1 = 0.0;
float	 BLUE_1 = 0.1;

float	  RED_2 = 0.0;
float	GREEN_2 = 0.3;
float	 BLUE_2 = 0.8;

float	  RED_3 = 1.0;
float	GREEN_3 = 1.0;
float	 BLUE_3 = 0.0;

float	  RED_4 = 0.8;
float	GREEN_4 = 0.3;
float	 BLUE_4 = 0.0;

float	  RED_5 = 0.1;
float	GREEN_5 = 0.0;
float	 BLUE_5 = 0.0;


void    free_vector();
void    free_matrix();
void    free_f3tensor();
void    free_int3tensor();
float   **matrix( long nrl, long nrh, long ncl, long nch );
float   ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh );
unsigned short int      ***int3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh );
void    nrerror( char error_text[]);                                                               
float   *vector();

void	do_plot_xy();
void 	two_columns();
void 	draw_density();
void 	do_draw_density();
void 	three_columns();
void	do_plot_two();
void	do_plot_contours();
void	contours();
float 	bicubic();
void	getcol( float val, int *r, int *g, int *b );
void	kzfilter();
int	HAVE_DONE_KZ = NO;
int	DRAW_KZ = NO;
int	halfwidth;


float	*x;
float	*y;
float	*z;
float	*kz;
float	*x_inter;
float	*y_inter;
float	*z_inter;
float	maxx, maxy;
float	minx, miny;
float	maxz, minz;
float	FIRST = 0;
float	dx, dy;
int	N, i;
int     DRAW_DENSITY = NO;
int	DRAW_HISTOGRAM = NO;
int	SAVE_HIST = NO;
int	CUMM_HIST = NO;
int	HAVE_LIMITS = NO;
long    winid;
char	title[300];
int	columns=0;
int	lines=0;
float	min, max, mean, rmsd;
float	c_step;
int	dots = NO;
int	filled_dots = NO;
float	**mat;
float	**mat_inter;
float	**inter;
float	**inter_inter;
int	VERBOSE = NO;
int	DRAW_NEG = NO;
int	IN_ZOOM = NO;
int	COLOR = NO;
int	DRAW_CONT = NO;
int	HAVE_COL = 0;
int	COL1 = 1;
int	COL2 = 2;
int	COL3 = 3;
int	SUBMATRIX = NO;
int	newc, newl;
int	oric, oril;

int	BICUB = YES;
int	STEP  = 500000;
int	MAXN  = 500000;
int	MAXML = 2000;
int	STEPM = 1000;

int	LOG = NO;
int	LOGLOG = NO;
int	PLOT_LABELS = 0;
char	label[100];
char	label2[100];

int	NOW_PLAYING = -1;
int	STEP_PLAYING = 1;
int	TIMING = NO;
time_t  start_timer;
time_t  end_timer;  


float	MAX_FROM_CML;
float	MIN_FROM_CML;
int	HAVE_MIN_MAX = NO;
int	HalfWidthOffset = 0;

int	COLOR_DOTS = NO;
float	hist_step;
int	PLOT_HIST = 1;
int	ERROR_BARS = NO;

int	FIRST_REDRAW = YES;
int	FINE_DENSITY_GRID = NO;
int	have_plot_labels = 0;
int	AUTOSCALE = NO;
float	scalesecond;
float	transsecond;

float	welford_mean;
float	welford_variance;
int	KZ_ROUNDS = 0;
float	KZ_ROUNDS_DS = 0;

void 	DFT(int m,float *y1);
int	HAVE_DFT = 0;
float	*data_DFT1;
float	*data_DFT2;

int	PLOT_LIGHTGRID = 0;
int     DATA_SET = -1;
int     HIST_EXACT = 0;
int     LARGE_LABELS = 0;

int     RUNAWK = NO;
char    tempfilename[1000] = "/tmp/plot_XXXXXX";
int     HAVE_TEMPFILE = NO;

int     SECTIONS = 1;
int     SECTION = 0;
int	mult = 0;
int	ori_col, ori_lin;

int     MINX_SIZE = 700;
int     MINY_SIZE = 300;
int     MINS_SIZE = 700;




int main(argc,argv)
int  	argc;
char	*argv[];
{
	short 	data;
	long	dev;
	char	line[500000];	
	float	junk;
	int	i, k;


        if ( getenv("PLOT_LABELS") )
          {
            PLOT_LABELS = 0;
            have_plot_labels = 0;
          }
        else
          {
            PLOT_LABELS = 1;
            have_plot_labels = 1;
          }
        
          
        if ( getenv("PLOT_DOTS") )
            dots = 1;
            
        if ( getenv("PLOT_HIST_EXACT") )
            HIST_EXACT = 1;

        if ( getenv("LARGE_LABELS") )
            {
                LARGE_LABELS = 1;
                MINX_SIZE = 1400;
                MINY_SIZE = 600;
                MINS_SIZE = 900;
            }
        
        if ( getenv("PLOT_LIGHTGRID") )
            PLOT_LIGHTGRID = 0;
        else
            PLOT_LIGHTGRID = 1;


	if ( argc >= 1 )
	{
	for ( i=1 ; i < argc ; i++ )
	{
		if ( strncasecmp( argv[i], "-V", 2 ) == 0 && strlen( argv[i] ) == 2 )
		{
			printf("This is plot v.%s\n", VERSION );
                        printf("\033[32m\033[1mDocs at : http://utopia.duth.gr/glykos/plot/\033[0m\n");
                        myexit( 0 );
		}
	}
	}

	if ( argc >= 1 )
	{
	for ( i=1 ; i < argc ; i++ )
	{
		if ( strncasecmp( argv[i], "-A", 2 ) == 0 && strlen( argv[i] ) == 2 )
		{
			AUTOSCALE = YES;
			for ( k=i+1 ; k < argc ; k++ )
				{
				strcpy( argv[k-1], argv[k] );
				}
			argc--;
		}
	}
	}


	if ( argc >= 1 )
	{
	for ( i=1 ; i < argc ; i++ )
	{
		if ( strncasecmp( argv[i], "-S", 2 ) == 0 && strlen( argv[i] ) == 2 )
		{
                        if ( sscanf( argv[i+1], "%d", &SECTIONS ) != 1 || SECTIONS < 1 )
                            {
			        printf("\033[31m\033[1mExpected a positive integer after the -s flag.\033[0m\n");
			        myexit(1);
                            }
			for ( k=i+1 ; k < argc ; k++ )
				{
				strcpy( argv[k-1], argv[k] );
				}
                        argc--;
			for ( k=i+1 ; k < argc ; k++ )
				{
				strcpy( argv[k-1], argv[k] );
				}
			argc--;
		}
	}
	}


	if ( argc >= 1 )
	{
	for ( i=1 ; i < argc ; i++ )
	{
		if ( strncasecmp( argv[i], "-F", 2 ) == 0 && strlen( argv[i] ) == 2 )
		{
			FINE_DENSITY_GRID = YES;
			for ( k=i+1 ; k < argc ; k++ )
				{
				strcpy( argv[k-1], argv[k] );
				}
			argc--;
		}
	}
	}

	
        if ( argc >= 1 )
	{
	for ( i=1 ; i < argc ; i++ )
	{
		if ( strncasecmp( argv[i], "-D", 2 ) == 0 && strlen( argv[i] ) > 2 )
		{
                        if ( isdigit( argv[i][2] ) != 0 )
                            {
                                DATA_SET = argv[i][2] - '0';
                            }
                        else
                            {
                                DATA_SET = toupper( argv[i][2] ) - 'A' +10;
                            }

                        if ( DATA_SET < 1 || DATA_SET > 35 )
                            DATA_SET = 1;

                        
			for ( k=i+1 ; k < argc ; k++ )
				{
				strcpy( argv[k-1], argv[k] );
				}
			argc--;
		}
	}
	}



	
	                      /* min - max */
	if ( argc == 5 )
	{
	    if ( strncasecmp( argv[1], "-R", 2 ) == 0 )
	      {
	        if ( sscanf( argv[2], "%f", &MIN_FROM_CML ) == 1 && sscanf( argv[3], "%f", &MAX_FROM_CML ) == 1 )
	        {
	          strcpy( argv[1], argv[4] );
	          argc = 2;
	          HAVE_MIN_MAX = YES;
	        }
	      }
	}

	if ( argc == 4 )
	{
	    if ( strncasecmp( argv[1], "-R", 2 ) == 0 )
	      {
	        if ( sscanf( argv[2], "%f", &MIN_FROM_CML ) == 1 && sscanf( argv[3], "%f", &MAX_FROM_CML ) == 1 )
	        {
	          argc = 1;
	          HAVE_MIN_MAX = YES;
	        }
	      }
	}
	
	
	if ( argc > 2 )
		{
			printf("\033[31m\033[1mToo many command line flags.\nForgotten input redirection maybe ?\033[0m\n");
			printf("\033[32m\033[1mDocs at : http://utopia.duth.gr/glykos/plot/\033[0m\n");
			myexit(1);
		}


	if ( argc == 2 )
		{
		        if ( strncasecmp( argv[1], "-LL", 3 ) == 0 )
		          {
		            LOGLOG = YES;
		            if ( strlen( argv[1] ) == 3 )
		            	argc = 1;
			    else
			    {
		            strcpy( line, &argv[1][3] );
		            strcpy( &argv[1][1], line );
		            }
		          }
		        else if ( strncasecmp( argv[1], "-L", 2 ) == 0 )
		          {
		            LOG = YES;
		            if ( strlen( argv[1] ) == 2 )
		            	argc = 1;
			    else
			    {
		            strcpy( line, &argv[1][2] );
		            strcpy( &argv[1][1], line );
		            }
		          }
		        
		        if ( argc == 1 )
		        	;
			else if ( strncasecmp( argv[1], "-HS", 3 ) == 0 && strlen( argv[1] ) == 3 )
				{
					DRAW_HISTOGRAM = YES;
					SAVE_HIST = YES;
					argc = 1;
				}
			else
			if ( strncasecmp( argv[1], "-HC", 3 ) == 0 && strlen( argv[1] ) == 3 )
				{
					DRAW_HISTOGRAM = YES;
					CUMM_HIST = YES;
					argc = 1;
				}
			else
			if ( strncasecmp( argv[1], "-HCS", 4 ) == 0 && strlen( argv[1] ) == 4 )
				{
					DRAW_HISTOGRAM = YES;
					CUMM_HIST = YES;
                                        SAVE_HIST = YES;
					argc = 1;
				}
			else
			if ( strncasecmp( argv[1], "-H", 2 ) == 0 && strlen( argv[1] ) == 2 )
				{
					DRAW_HISTOGRAM = YES;
					argc = 1;
				}
			else
			if ( strncasecmp( argv[1], "-F", 2 ) == 0 && strlen( argv[1] ) == 2 )
				{
					FINE_DENSITY_GRID = YES;
					argc = 1;
				}
			else
			if ( strncasecmp( argv[1], "-XYDY", 5 ) == 0 && strlen( argv[1] ) == 5 )
				{
					ERROR_BARS = YES;
					argc = 1;
				}
			else
			if ( strncasecmp( argv[1], "-CV", 3 ) == 0 )
				{
					if ( strncasecmp( argv[1], "-CVS", 4 ) == 0 )
					{
					SUBMATRIX = YES;
					if ( sscanf( &argv[1][4], "%d,%d,%d,%d", &newc, &newl, &oric, &oril) != 4 )
						{
						printf("\033[31m\033[1mWrong arguments for submatrix.\033[0m\n");
						printf("\033[32m\033[1mDocs at : http://utopia.duth.gr/glykos/plot/\033[0m\n");
						myexit(1);
						}
					}
						
					VERBOSE = YES;
                                        contours();
				}
			else
			if ( strncasecmp( argv[1], "-CC", 3 ) == 0 )
				{
					if ( strncasecmp( argv[1], "-CCS", 4 ) == 0 )
					{
					SUBMATRIX = YES;
					if ( sscanf( &argv[1][4], "%d,%d,%d,%d", &newc, &newl, &oric, &oril) != 4 )
						{
						printf("\033[31m\033[1mWrong arguments for submatrix.\033[0m\n");
						printf("\033[32m\033[1mDocs at : http://utopia.duth.gr/glykos/plot/\033[0m\n");
						myexit(1);
						}
					}
					                                                
					COLOR = YES;
                                        contours();
				}
			else
			if ( strncasecmp( argv[1], "-C", 2 ) == 0 )
				{
					if ( strncasecmp( argv[1], "-CS", 3 ) == 0 )
					{
					SUBMATRIX = YES;
					if ( sscanf( &argv[1][3], "%d,%d,%d,%d", &newc, &newl, &oric, &oril) != 4 )
						{
						printf("\033[31m\033[1mWrong arguments for submatrix.\033[0m\n");
						printf("\033[32m\033[1mDocs at : http://utopia.duth.gr/glykos/plot/\033[0m\n");
						myexit(1);
						}
					}
					                                                
                                        contours();
				}
			else
			if ( 	strncasecmp( argv[1], "-K", 2 ) == 0 )		/* define columns from the command line */
				{
                                        if ( strncmp( argv[1], "-K", 2 ) == 0 )
                                            RUNAWK = YES;

					if ( strlen( argv[1] ) <= 2 || strlen( argv[1] ) >= 6 )
					{
					printf("\033[31m\033[1mColumn definition is wrong.\033[0m\n");
					printf("\033[32m\033[1mDocs at : http://utopia.duth.gr/glykos/plot/\033[0m\n");
					myexit(1);
					}

										/* We sure have one column */
					HAVE_COL = 1;

					if ( isdigit( argv[1][2] ) != 0 )
						{
							COL1 = argv[1][2] - '0';
						}
					else
						{
							COL1 = toupper( argv[1][2] ) - 'A' +10;
						}

					
					if ( strlen( argv[1] ) > 3 )		/* Two columns */
						{
							HAVE_COL = 2;
							
							if ( isdigit( argv[1][3] ) != 0 )
								{
									COL2 = argv[1][3] - '0';
								}
							else
								{
									COL2 = toupper( argv[1][3] ) - 'A' +10;
								}

						}

					if ( strlen( argv[1] ) == 5 )		/* Three columns */
						{
							HAVE_COL = 3;
							
							if ( isdigit( argv[1][4] ) != 0 )
								{
									COL3 = argv[1][4] - '0';
								}
							else
								{
									COL3 = toupper( argv[1][4] ) - 'A' +10;
								}

						}


					if ( COL1 < 1 || COL2 < 1 || COL3 < 1 )
					{
					printf("\033[31m\033[1mColumn definition is wrong.\033[0m\n");
					printf("\033[32m\033[1mDocs at : http://utopia.duth.gr/glykos/plot/\033[0m\n");
					myexit(1);
					}
					
					argc = 1;
				}
			else
			if ( 	strncasecmp( argv[1], "-HK", 3 ) == 0 )		/* hist + define columns from the command line */
				{
					if ( strlen( argv[1] ) <= 3 || strlen( argv[1] ) >= 6 )
					{
					printf("\033[31m\033[1mColumn definition is wrong.\033[0m\n");
					printf("\033[32m\033[1mDocs at : http://utopia.duth.gr/glykos/plot/\033[0m\n");
					myexit(1);
					}

										/* We sure have one column */
					HAVE_COL = 1;

					if ( isdigit( argv[1][3] ) != 0 )
						{
							COL1 = argv[1][3] - '0';
						}
					else
						{
							COL1 = toupper( argv[1][3] ) - 'A' +10;
						}

					
					if ( strlen( argv[1] ) == 5 )		/* Two columns */
						{
							HAVE_COL = 2;
							
							if ( isdigit( argv[1][4] ) != 0 )
								{
									COL2 = argv[1][4] - '0';
								}
							else
								{
									COL2 = toupper( argv[1][4] ) - 'A' +10;
								}

						}

					if ( strlen( argv[1] ) > 5 )		/* Three columns */
						{
						printf("\033[31m\033[1mToo many columns for histogram.\033[0m\n");
						printf("\033[32m\033[1mDocs at : http://utopia.duth.gr/glykos/plot/\033[0m\n");
						myexit(1);
						}


					if ( COL1 < 1 || COL2 < 1 )
					{
					printf("\033[31m\033[1mColumn definition is wrong.\033[0m\n");
					printf("\033[32m\033[1mDocs at : http://utopia.duth.gr/glykos/plot/\033[0m\n");
					myexit(1);
					}
					
					DRAW_HISTOGRAM = YES;
					argc = 1;
				}
			else
				{

					printf("\033[31m\033[1mCommand line flag does not make sense.\033[0m\n");
					printf("\033[31m\033[1mForgotten input redirection maybe ?\033[0m\n");
					printf("\033[32m\033[1mDocs at : http://utopia.duth.gr/glykos/plot/\033[0m\n");
					myexit(1);
				}
		}

	
	if ( argc != 1 )
		{
			printf("\033[31m\033[1mThis is a unix filter. Use redirection or pipes.\033[0m\n");
			printf("\033[32m\033[1mDocs at : http://utopia.duth.gr/glykos/plot/\033[0m\n");
			myexit(1);
		}
		


        if ( RUNAWK == NO )
        {
            i = 0;
            while( fgets( line, 499999, stdin ) != NULL )
            {
            if ( (columns = sscanf( line, "%f %f %f %f", &junk, &junk, &junk, &junk)) >= 1 )
                    break;
            
            i++;
            if ( i < 2 )
            printf("\033[37m\033[1mCaution: header line skipped:\033[0m %s", line );	
            }	
            if ( i > 1 )
                printf("\033[37m\033[1m     ... followed by %d more lines.\033[0m\n", i );
        }

        /* 
         * Attempt to deal with lazy users that have non-numerical columns
         *
         */

        if ( (HAVE_COL > 0 && columns < HAVE_COL && RUNAWK == NO) || (HAVE_COL > 0 && RUNAWK == YES) )
        {
            char    scall[1000];
            int     sysret;

            if ( RUNAWK == NO )
            {
                printf("\033[32m\033[1mAre you feeding plot non-numerical data ?!? Beware ...\033[0m\n");
                rewind( stdin );
            }
            else
            {
                printf("\033[32m\033[1mWill attempt to use awk to select columns. This may not go as planned ...\033[0m\n");
            }

            close( mkstemp( tempfilename ));
            HAVE_TEMPFILE = YES;

            if ( HAVE_COL == 1 )
            {
                sprintf( scall, "awk '{print $%d}' /dev/stdin > %s", COL1, tempfilename ); 
                COL1 = 1;
            }
            if ( HAVE_COL == 2 )
            {
                sprintf( scall, "awk '{print $%d,$%d}' /dev/stdin > %s", COL1, COL2, tempfilename ); 
                COL1 = 1;
                COL2 = 2;
            }
            if ( HAVE_COL == 3 )
            {
                sprintf( scall, "awk '{print $%d,$%d,$%d}' /dev/stdin > %s", COL1, COL2, COL3, tempfilename ); 
                COL1 = 1;
                COL2 = 2;
                COL2 = 3;
            }
            sysret = system( scall );
            if ( sysret < 0 )
            {
		printf("\033[31m\033[1mFailed to execute awk. This is all that is known.\033[0m\n");
    		myexit(1);
            }
            stdin = freopen( tempfilename , "r", stdin);

            while( fgets( line, 499999, stdin ) != NULL )
            {
            if ( (columns = sscanf( line, "%f %f %f %f", &junk, &junk, &junk, &junk)) >= 1 )
                    break;
            
            printf("\033[37m\033[1mCaution: header line skipped:\033[0m %s", line );	
            }	

        }



	if ( columns > 3 )
		{
			columns = -2;
		}


        x = vector( 0, MAXN );
        y = vector( 0, MAXN );
        z = vector( 0, MAXN );


	/* One column, do a simple x-y plot */
	if ( columns == 1 || HAVE_COL == 1  )
	{

	columns = 1;

	if ( HAVE_COL != 0 )
	{
	char 	*p;
	int	tot_col;
	float	*mat;
	int	k;
	float	val;


	
	/* Number of columns */
	p = &line[0];
	tot_col = 0;
	while ( sscanf( p, "%f%n", &junk, &i) == 1 )
		{
			tot_col++;
			p += i;
		}
	
	mat = vector( 0, tot_col);

	if ( COL1 > tot_col )
		{
			printf("\033[31m\033[1mNo such column (%d). Max is %d.\033[0m\n", COL1, tot_col);
			printf("\033[32m\033[1mDocs at : http://utopia.duth.gr/glykos/plot/\033[0m\n");
			myexit(1);		
		}

	/* Grab first number */
	p = &line[0];
	x[0] = 1;
	k = 0;
	while ( sscanf( p, "%f%n", &y[0], &i) == 1 )
		{
			if ( COL1 == (k+1))
				break;
			k++;
			p += i;
		}



	/* Read the rest */
	N = 1;
	while( 1 )
	{
		for ( i=0 ; i < tot_col ; i++ )
			{
				if ( scanf("%f", &val) != 1 )
					break;
                                if ( isnan(val) || isinf(val) )
                                  {
                                    printf("\033[31m\033[1mData appear to contain one or more 'NaN' and/or 'Inf'. Goodbye.\033[0m\n");
                                    myexit(1);
                                  }
                                
				mat[i] = val;
			}
		
		if ( i == 0 )
			break;
		else if ( i != tot_col )
			{
				printf("\033[31m\033[1mNumber of columns in matrix not constant ? Abort.\033[0m\n");
				printf("\033[32m\033[1mDocs at : http://utopia.duth.gr/glykos/plot/\033[0m\n");
				myexit(1);
			}

		x[N] = N+1;
		y[N] = mat[ COL1-1 ];
		N++;

		if ( N == MAXN-1 )
			{
			    x_inter = vector( 0, MAXN+STEP );
			    y_inter = vector( 0, MAXN+STEP );
			    z_inter = vector( 0, MAXN+STEP );
                            memcpy( x_inter, x, (MAXN+1)*sizeof(float) );
                            memcpy( y_inter, y, (MAXN+1)*sizeof(float) );
                            memcpy( z_inter, z, (MAXN+1)*sizeof(float) );
                            free_vector( x, 0, MAXN );
                            free_vector( y, 0, MAXN );
                            free_vector( z, 0, MAXN );
                            x = x_inter;
                            y = y_inter;
                            z = z_inter;
                            MAXN += STEP;
			}
			

	}



	}
	else
	{

	sscanf( line, "%f", &y[0] );
	x[0] = 1;
	N = 1;
	while( scanf( "%f", &y[N] ) == 1 )
		{
			x[N] = N+1;
			N++;
			if ( N == MAXN-1 )
				{
                                    x_inter = vector( 0, MAXN+STEP );
			            y_inter = vector( 0, MAXN+STEP );
			            z_inter = vector( 0, MAXN+STEP );
                                    memcpy( x_inter, x, (MAXN+1)*sizeof(float) );
                                    memcpy( y_inter, y, (MAXN+1)*sizeof(float) );
                                    memcpy( z_inter, z, (MAXN+1)*sizeof(float) );
                                    free_vector( x, 0, MAXN );
                                    free_vector( y, 0, MAXN );
                                    free_vector( z, 0, MAXN );
                                    x = x_inter;
                                    y = y_inter;
                                    z = z_inter;
                                    MAXN += STEP;
				}
		}

	}


	if ( LOGLOG == YES )
	  {
	    for ( i=0 ; i < N ; i++ )
	      {
	        if ( y[i] <= 0 || x[i] <= 0 )
			{
				printf("\033[31m\033[1mLogarithm requested but negative values present ? Abort.\033[0m\n");
				myexit(1);
			}
	        y[i] = log( y[i] );
	        x[i] = log( x[i] );
	      }
          }
	if ( LOG == YES )
	  {
	    for ( i=0 ; i < N ; i++ )
	      {
	        if ( y[i] <= 0 )
			{
				printf("\033[31m\033[1mLogarithm requested but negative values present ? Abort.\033[0m\n");
				myexit(1);
			}
	        y[i] = log( y[i] );
	      }
          }

        if ( N > 20000 && DRAW_HISTOGRAM == NO )
        {   
            dots = 1;
        }

	if ( DRAW_HISTOGRAM == NO )
		two_columns();
	else
	{
		int	h_N;
		double	max, min;
		double	step;


		
		max = y[0];
		min = y[0];
		
		for ( i=0 ; i < N ; i++ )
			{
				if ( y[i] > max )
					max = y[i];
				if ( y[i] < min )
					min = y[i];
			}


		/* Freedman-Diaconis rule */
		sort( N, y-1 );
		h_N = (int)( (max-min) / (2*(y[ (int)(3.0*N/4.0 + 0.5)]-y[ (int)(N/4.0 +0.5) ]) / pow( N, 1.0l/3.0l)) + 0.50 );

		if ( h_N < 2 )
		{
			h_N = (int)( 2.0 * sqrt(N) + 0.50);
			if ( h_N < 2 )
			{
				printf("Not enough data for drawing a histogram. Goodbye.\n\n");
				myexit(1);
			}
			fprintf(stderr, "\033[31m\033[1mFreedman-Diaconis rule failed. Be skeptical.\033[0m\n");			
		}

		if ( h_N < 2 )
			{
				printf("Not enough data for drawing a histogram. Goodbye.\n\n");
				myexit(1);
			}
			
		if ( h_N > 500 )
			{
			fprintf(stderr, "\033[31m\033[1mToo many bins from Freedman-Diaconis rule. Set to 500 bins.\033[0m\n");
			h_N = 500;
			}
		
		step = ( max - min ) / h_N ;

		if ( HAVE_MIN_MAX == YES )
	        {
	        step = MIN_FROM_CML;
	        if ( step <= 0 )
	        {
	            printf("\033[31m\033[1mThe width chosen for the histogram bins must be greater than zero.\033[0m\n");
	            myexit(1);
	        }
	        h_N  = (int)((max-min) / step + 0.50 );
	        if ( h_N > 999 )
	          {
	            printf("\033[31m\033[1mThe width chosen for the histogram bins is too small. Increase.\033[0m\n");
	            myexit(1);
	          }
	        step = (max-min) / h_N;
	        }


                if ( HIST_EXACT == 0 )
                {
	        h_N++;
	        min = min - step / 2 ;
	        max = max + step / 2 ;
                }

		if ( step == 0 )
			{
				printf("\033[31m\033[1mConstant y values ? Goodbye.\033[0m\n");
				myexit(1);
			}

                x_inter = vector( 0, MAXN+STEP );
		y_inter = vector( 0, MAXN+STEP );
		z_inter = vector( 0, MAXN+STEP );
                memcpy( x_inter, x, (MAXN+1)*sizeof(float) );
                memcpy( y_inter, y, (MAXN+1)*sizeof(float) );
                memcpy( z_inter, z, (MAXN+1)*sizeof(float) );
                free_vector( x, 0, MAXN );
                free_vector( y, 0, MAXN );
                free_vector( z, 0, MAXN );
                x = x_inter;
                y = y_inter;
                z = z_inter;
                MAXN += STEP;

		
		if ( N > MAXN - 1000 )
			N = MAXN - 1000;
		
		for ( i=0 ; i < 1000 ; i++ )
			{
				x[ MAXN - 1000 + i] = 0.0;
				y[ MAXN - 1000 + i] = 0.0;
			}
		
		for ( i=0 ; i < h_N ; i++ )
			x[ MAXN - 1000 + i] = min + step*( i + 0.50 );
		
		for ( i=0 ; i < N ; i++ )
			y[ MAXN - 1000 + (int)((y[i] - min) / step) ]++;

                if ( HIST_EXACT == 1 )
                    y[ MAXN - 1000 + h_N - 1] += y[ MAXN - 1000 + h_N ];

		for ( i=0 ; i < h_N ; i++ )
			{
				x[i] = x[ MAXN - 1000 + i];
				y[i] = y[ MAXN - 1000 + i];
			}
		N = h_N;

		if ( CUMM_HIST == YES )
		  {
		    float prev = 0.0;
		    
		    for ( i=0 ; i < h_N ; i++ )
			{
				y[i] = prev + y[i];
				prev = y[i];
			}
                    
                    if ( AUTOSCALE == YES )
                    {
                        float   maxcum;

                        maxcum = y[i-1];
                        for ( i=0 ; i < h_N ; i++ )
                            y[i] /= maxcum;

                    }

		  }

		if ( SAVE_HIST == YES)
		{
		  FILE *histf;
		  
		  histf = fopen( "plot.histogram", "w" );
		  if ( histf != NULL )
		  {
		    for ( i=0 ; i < h_N ; i++ )
		      fprintf( histf, " % 12.10e % 12.10e\n", x[i], y[i] );
		    
		    fclose( histf );
		  }
		}
		
		hist_step = step;
		two_columns();
			
	}
	
	}



	
	
	/* Two columns, do a simple x-y plot */
	if ( columns == 2 || HAVE_COL == 2 )
	{


	columns = 2;

	if ( HAVE_COL != 0 )
	{
	char 	*p;
	int	tot_col;
	float	*mat;
	int	k;
	float	val;
	
	/* Number of columns */
	p = &line[0];
	tot_col = 0;
	while ( sscanf( p, "%f%n", &junk, &i) == 1 )
		{
			tot_col++;
			p += i;
		}
	
	mat = vector( 0, tot_col);


	if ( COL1 > tot_col || COL2 > tot_col )
		{
			printf("\033[31m\033[1mNo such column. Max is %d.\033[0m\n", tot_col);
			myexit(1);		
		}


	/* Grab first number */
	p = &line[0];
	k = 0;
	while ( sscanf( p, "%f%n", &mat[k], &i) == 1 )
		{
			k++;
			p += i;
		}
	x[0] = mat[ COL1 - 1 ];
	y[0] = mat[ COL2 - 1 ];



	/* Read the rest */
	N = 1;
	while( 1 )
	{
		for ( i=0 ; i < tot_col ; i++ )
			{
				if ( scanf("%f", &val) != 1 )
					break;
                                if ( isnan(val) || isinf(val) )
                                  {
                                    printf("\033[31m\033[1mData appear to contain one or more 'NaN' and/or 'Inf'. Goodbye.\033[0m\n");
                                    myexit(1);
                                  }
				mat[i] = val;
			}
		
		if ( i == 0 )
			break;
		else if ( i != tot_col )
			{
				printf("\033[31m\033[1mNumber of columns in matrix not constant ? Abort.\033[0m\n");
				myexit(1);
			}

		x[N] = mat[ COL1-1 ];
		y[N] = mat[ COL2-1 ];
		N++;

		if ( N == MAXN-1 )
			{
                                    x_inter = vector( 0, MAXN+STEP );
			            y_inter = vector( 0, MAXN+STEP );
			            z_inter = vector( 0, MAXN+STEP );
                                    memcpy( x_inter, x, (MAXN+1)*sizeof(float) );
                                    memcpy( y_inter, y, (MAXN+1)*sizeof(float) );
                                    memcpy( z_inter, z, (MAXN+1)*sizeof(float) );
                                    free_vector( x, 0, MAXN );
                                    free_vector( y, 0, MAXN );
                                    free_vector( z, 0, MAXN );
                                    x = x_inter;
                                    y = y_inter;
                                    z = z_inter;
                                    MAXN += STEP;
			}
			

	}



	}
	else
	{



	sscanf( line, "%f %f", &x[0], &y[0] );
	N = 1;
	while( scanf( "%f %f", &x[N], &y[N] ) == 2 )
		{
			if ( isnan(x[N]) || isinf(x[N]) || isnan(y[N]) || isinf(y[N]) )
				{
                                printf("\033[31m\033[1mData appear to contain one or more 'NaN' and/or 'Inf'. Goodbye.\033[0m\n");
                                myexit(1);
                                }
			N++;
			if ( N == MAXN-1 )
				{
                                    x_inter = vector( 0, MAXN+STEP );
			            y_inter = vector( 0, MAXN+STEP );
			            z_inter = vector( 0, MAXN+STEP );
                                    memcpy( x_inter, x, (MAXN+1)*sizeof(float) );
                                    memcpy( y_inter, y, (MAXN+1)*sizeof(float) );
                                    memcpy( z_inter, z, (MAXN+1)*sizeof(float) );
                                    free_vector( x, 0, MAXN );
                                    free_vector( y, 0, MAXN );
                                    free_vector( z, 0, MAXN );
                                    x = x_inter;
                                    y = y_inter;
                                    z = z_inter;
                                    MAXN += STEP;
				}
		}
	}




        if ( N > 20000 && DRAW_HISTOGRAM == NO )
        {   
            dots = 1;
        }


	if ( LOGLOG == YES )
	  {
	    for ( i=0 ; i < N ; i++ )
	      {
	        if ( y[i] <= 0 || x[i] <= 0 )
			{
				printf("\033[31m\033[1mLogarithm requested but negative values present ? Abort.\033[0m\n");
				myexit(1);
			}
	        y[i] = log( y[i] );
	        x[i] = log( x[i] );
	      }
          }
	if ( LOG == YES )
	  {
	    for ( i=0 ; i < N ; i++ )
	      {
	        if ( y[i] <= 0 )
			{
				printf("\033[31m\033[1mLogarithm requested but negative values present ? Abort.\033[0m\n");
				myexit(1);
			}
	        y[i] = log( y[i] );
	      }
          }



	if ( DRAW_HISTOGRAM == NO )
		two_columns();
	else
	{
		int	h_N;
		double	max, min;
		double	step;


		max = y[0];
		min = y[0];
		
		for ( i=0 ; i < N ; i++ )
			{
				if ( y[i] > max )
					max = y[i];
				if ( y[i] < min )
					min = y[i];
			}

                

		/* Freedman-Diaconis rule */
		sort( N, y-1 );
		h_N = (int)( (max-min) / (2*(y[ (int)(3.0*N/4.0 + 0.5)]-y[ (int)(N/4.0 +0.5) ]) / pow( N, 1.0l/3.0l)) + 0.50 );

		if ( h_N < 2 )
		{
			h_N = (int)( 2.0 * sqrt(N) + 0.50);
			if ( h_N < 2 )
			{
				printf("Not enough data for drawing a histogram. Goodbye.\n\n");
				myexit(1);
			}
			fprintf(stderr, "\033[31m\033[1mFreedman-Diaconis rule failed. Be skeptical.\033[0m\n");			
		}

		if ( h_N < 2 )
			{
				printf("Not enough data for drawing a histogram. Goodbye.\n\n");
				myexit(1);
			}
			
		if ( h_N > 500 )
			{
			fprintf(stderr, "\033[31m\033[1mToo many bins from Freedman-Diaconis rule. Set to 500 bins.\033[0m\n");
			h_N = 500;
			}
		
		step = ( max - min ) / h_N ;

		if ( HAVE_MIN_MAX == YES )
	        {
	        step = MIN_FROM_CML;
	        h_N  = (int)((max-min) / step + 0.50 );
	        if ( h_N > 999 )
	          {
	            printf("\033[31m\033[1mThe width chosen for the histogram bins is too small. Increase.\033[0m\n");
	            myexit(1);
	          }
	        step = (max-min) / h_N;
	        }

                if ( HIST_EXACT == 0 )
                {
	        h_N++;
	        min = min - step / 2 ;
	        max = max + step / 2 ;
                }

		if ( step == 0 )
			{
				printf("\033[31m\033[1mConstant y values ? Goodbye.\033[0m\n");
				myexit(1);
			}

                x_inter = vector( 0, MAXN+STEP );
		y_inter = vector( 0, MAXN+STEP );
		z_inter = vector( 0, MAXN+STEP );
                memcpy( x_inter, x, (MAXN+1)*sizeof(float) );
                memcpy( y_inter, y, (MAXN+1)*sizeof(float) );
                memcpy( z_inter, z, (MAXN+1)*sizeof(float) );
                free_vector( x, 0, MAXN );
                free_vector( y, 0, MAXN );
                free_vector( z, 0, MAXN );
                x = x_inter;
                y = y_inter;
                z = z_inter;
                MAXN += STEP;

		
		if ( N > MAXN - 1000 )
			N = MAXN - 1000;
		
		for ( i=0 ; i < 1000 ; i++ )
			{
				x[ MAXN - 1000 + i] = 0.0;
				y[ MAXN - 1000 + i] = 0.0;
			}
		
		for ( i=0 ; i < h_N ; i++ )
			x[ MAXN - 1000 + i] = min + step*( i + 0.50 );
		
		for ( i=0 ; i < N ; i++ )
			y[ MAXN - 1000 + (int)((y[i] - min) / step) ]++;
			
                if ( HIST_EXACT == 1 )
                    y[ MAXN - 1000 + h_N - 1] += y[ MAXN - 1000 + h_N ];

		for ( i=0 ; i < h_N ; i++ )
			{
				x[i] = x[ MAXN - 1000 + i];
				y[i] = y[ MAXN - 1000 + i];
			}

		N = h_N;

		if ( CUMM_HIST == YES )
		  {
		    float prev = 0.0;
		    
		    for ( i=0 ; i < h_N ; i++ )
			{
				y[i] = prev + y[i];
				prev = y[i];
			}
                    
                    if ( AUTOSCALE == YES )
                    {
                        float   maxcum;

                        maxcum = y[i-1];
                        for ( i=0 ; i < h_N ; i++ )
                            y[i] /= maxcum;

                    }
		  }

		if ( SAVE_HIST == YES)
		{
		  FILE *histf;
		  
		  histf = fopen( "plot.histogram", "w" );
		  if ( histf != NULL )
		  {
		    for ( i=0 ; i < h_N ; i++ )
		      fprintf( histf, " % 12.10e % 12.10e\n", x[i], y[i] );
		    
		    fclose( histf );
		  }
		}
		
		hist_step = step;
		two_columns();
			
	}

	}




	
	/* Too many columns, do a simple x-y plot of the first two columns */
	if ( columns == -2 && HAVE_COL == 0 )
	{
	sscanf( line, "%f %f", &x[0], &y[0] );
	N = 1;
	while ( fgets( line, 499999, stdin ) != NULL )
		{
			if ( sscanf( line, "%f %f", &x[N], &y[N] ) == 2 )
			{
                                
			if ( isnan(x[N]) || isinf(x[N]) || isnan(y[N]) || isinf(y[N]) )
				{
                                printf("\033[31m\033[1mData appear to contain one or more 'NaN' and/or 'Inf'. Goodbye.\033[0m\n");
                                myexit(1);
                                }

			N++;
			if ( N == MAXN-1 )
				{
                                    x_inter = vector( 0, MAXN+STEP );
			            y_inter = vector( 0, MAXN+STEP );
			            z_inter = vector( 0, MAXN+STEP );
                                    memcpy( x_inter, x, (MAXN+1)*sizeof(float) );
                                    memcpy( y_inter, y, (MAXN+1)*sizeof(float) );
                                    memcpy( z_inter, z, (MAXN+1)*sizeof(float) );
                                    free_vector( x, 0, MAXN );
                                    free_vector( y, 0, MAXN );
                                    free_vector( z, 0, MAXN );
                                    x = x_inter;
                                    y = y_inter;
                                    z = z_inter;
                                    MAXN += STEP;
				}
                        }
                        else
                        {
                          printf("\033[37m\033[1mCaution: line skipped:\033[0m %s", line );
                        }
		}


        if ( N > 20000 && DRAW_HISTOGRAM == NO )
        {   
            dots = 1;
        }

	if ( LOGLOG == YES )
	  {
	    for ( i=0 ; i < N ; i++ )
	      {
	        if ( y[i] <= 0 || x[i] <= 0 )
			{
				printf("\033[31m\033[1mLogarithm requested but negative values present ? Abort.\033[0m\n");
				myexit(1);
			}
	        y[i] = log( y[i] );
	        x[i] = log( x[i] );
	      }
          }
	if ( LOG == YES )
	  {
	    for ( i=0 ; i < N ; i++ )
	      {
	        if ( y[i] <= 0 )
			{
				printf("\033[31m\033[1mLogarithm requested but negative values present ? Abort.\033[0m\n");
				myexit(1);
			}
	        y[i] = log( y[i] );
	      }
          }



	if ( DRAW_HISTOGRAM == NO )
		two_columns();
	else
	{
		int	h_N;
		double	max, min;
		double	step;

		max = y[0];
		min = y[0];
		
		for ( i=0 ; i < N ; i++ )
			{
				if ( y[i] > max )
					max = y[i];
				if ( y[i] < min )
					min = y[i];
			}


		/* Freedman-Diaconis rule */
		sort( N, y-1 );
		h_N = (int)( (max-min) / (2*(y[ (int)(3.0*N/4.0 + 0.5)]-y[ (int)(N/4.0 +0.5) ]) / pow( N, 1.0l/3.0l)) + 0.50 );

		if ( h_N < 2 )
		{
			h_N = (int)( 2.0 * sqrt(N) + 0.50);
			if ( h_N < 2 )
			{
				printf("Not enough data for drawing a histogram. Goodbye.\n\n");
				myexit(1);
			}
			fprintf(stderr, "\033[31m\033[1mFreedman-Diaconis rule failed. Be skeptical.\033[0m\n");			
		}

		if ( h_N < 2 )
			{
				printf("Not enough data for drawing a histogram. Goodbye.\n\n");
				myexit(1);
			}
			
		if ( h_N > 500 )
			{
			fprintf(stderr, "\033[31m\033[1mToo many bins from Freedman-Diaconis rule. Set to 500 bins.\033[0m\n");
			h_N = 500;
			}
		
		step = ( max - min ) / h_N ;

		if ( HAVE_MIN_MAX == YES )
	        {
	        step = MIN_FROM_CML;
	        h_N  = (int)((max-min) / step + 0.50 );
	        if ( h_N > 999 )
	          {
	            printf("\033[31m\033[1mThe width chosen for the histogram bins is too small. Increase.\033[0m\n");
	            myexit(1);
	          }
	        step = (max-min) / h_N;
	        }

                if ( HIST_EXACT == 0 )
                {
	        h_N++;
	        min = min - step / 2 ;
	        max = max + step / 2 ;
                }

		if ( step == 0 )
			{
				printf("\033[31m\033[1mConstant y values ? Goodbye.\033[0m\n");
				myexit(1);
			}

                x_inter = vector( 0, MAXN+STEP );
		y_inter = vector( 0, MAXN+STEP );
		z_inter = vector( 0, MAXN+STEP );
                memcpy( x_inter, x, (MAXN+1)*sizeof(float) );
                memcpy( y_inter, y, (MAXN+1)*sizeof(float) );
                memcpy( z_inter, z, (MAXN+1)*sizeof(float) );
                free_vector( x, 0, MAXN );
                free_vector( y, 0, MAXN );
                free_vector( z, 0, MAXN );
                x = x_inter;
                y = y_inter;
                z = z_inter;
                MAXN += STEP;

		
		if ( N > MAXN - 1000 )
			N = MAXN - 1000;
		
		for ( i=0 ; i < 1000 ; i++ )
			{
				x[ MAXN - 1000 + i] = 0.0;
				y[ MAXN - 1000 + i] = 0.0;
			}
		
		for ( i=0 ; i < h_N ; i++ )
			x[ MAXN - 1000 + i] = min + step*( i + 0.50 );
		
		for ( i=0 ; i < N ; i++ )
			y[ MAXN - 1000 + (int)((y[i] - min) / step) ]++;
			
                if ( HIST_EXACT == 1 )
                    y[ MAXN - 1000 + h_N - 1] += y[ MAXN - 1000 + h_N ];

		for ( i=0 ; i < h_N ; i++ )
			{
				x[i] = x[ MAXN - 1000 + i];
				y[i] = y[ MAXN - 1000 + i];
			}

		N = h_N;

		if ( CUMM_HIST == YES )
		  {
		    float prev = 0.0;
		    
		    for ( i=0 ; i < h_N ; i++ )
			{
				y[i] = prev + y[i];
				prev = y[i];
			}
                    
                    if ( AUTOSCALE == YES )
                    {
                        float   maxcum;

                        maxcum = y[i-1];
                        for ( i=0 ; i < h_N ; i++ )
                            y[i] /= maxcum;

                    }
		  }

		if ( SAVE_HIST == YES)
		{
		  FILE *histf;
		  
		  histf = fopen( "plot.histogram", "w" );
		  if ( histf != NULL )
		  {
		    for ( i=0 ; i < h_N ; i++ )
		      fprintf( histf, " % 12.10e % 12.10e\n", x[i], y[i] );
		    
		    fclose( histf );
		  }
		}
            
		hist_step = step;
		two_columns();
			
	}

	}






	/* Three columns, overlay two x-y plots, assuming first column is x */
	if ( columns == 3 || HAVE_COL == 3 )
	{


	columns = 3;

	if ( HAVE_COL != 0 )
	{
	char 	*p;
	int	tot_col;
	float	*mat;
	int	k;
	float	val;
	
	/* Number of columns */
	p = &line[0];
	tot_col = 0;
	while ( sscanf( p, "%f%n", &junk, &i) == 1 )
		{
			tot_col++;
			p += i;
		}
	
        mat = vector( 0, tot_col);
        
	if ( COL1 > tot_col || COL2 > tot_col || COL3 > tot_col )
		{
			printf("\033[31m\033[1mNo such column. Max is %d.\033[0m\n", tot_col);
			myexit(1);		
		}


	/* Grab first number */
	p = &line[0];
	k = 0;
	while ( sscanf( p, "%f%n", &mat[k], &i) == 1 )
		{
			k++;
			p += i;
		}
	x[0] = mat[ COL1 - 1 ];
	y[0] = mat[ COL2 - 1 ];
	z[0] = mat[ COL3 - 1 ];



	/* Read the rest */
	N = 1;
	while( 1 )
	{
		for ( i=0 ; i < tot_col ; i++ )
			{
				if ( scanf("%f", &val) != 1 )
					break;
                                if ( isnan(val) || isinf(val) )
                                  {
                                    printf("\033[31m\033[1mData appear to contain one or more 'NaN' and/or 'Inf'. Goodbye.\033[0m\n");
                                    myexit(1);
                                  }
				mat[i] = val;
			}
		
		if ( i == 0 )
			break;
		else if ( i != tot_col )
			{
				printf("\033[31m\033[1mNumber of columns in matrix not constant ? Abort.\033[0m\n");
				myexit(1);
			}

		x[N] = mat[ COL1-1 ];
		y[N] = mat[ COL2-1 ];
		z[N] = mat[ COL3-1 ];
		N++;

		if ( N == MAXN-1 )
			{
                                    x_inter = vector( 0, MAXN+STEP );
			            y_inter = vector( 0, MAXN+STEP );
			            z_inter = vector( 0, MAXN+STEP );
                                    memcpy( x_inter, x, (MAXN+1)*sizeof(float) );
                                    memcpy( y_inter, y, (MAXN+1)*sizeof(float) );
                                    memcpy( z_inter, z, (MAXN+1)*sizeof(float) );
                                    free_vector( x, 0, MAXN );
                                    free_vector( y, 0, MAXN );
                                    free_vector( z, 0, MAXN );
                                    x = x_inter;
                                    y = y_inter;
                                    z = z_inter;
                                    MAXN += STEP;
			}
			

	}



	}
	else
	{
	
	
	sscanf( line, "%f %f %f", &x[0], &y[0], &z[0] );
	N = 1;
	while( scanf( "%f %f %f", &x[N], &y[N], &z[N] ) == 3 )
		{
			N++;
			if ( N == MAXN-1 )
				{
                                    x_inter = vector( 0, MAXN+STEP );
			            y_inter = vector( 0, MAXN+STEP );
			            z_inter = vector( 0, MAXN+STEP );
                                    memcpy( x_inter, x, (MAXN+1)*sizeof(float) );
                                    memcpy( y_inter, y, (MAXN+1)*sizeof(float) );
                                    memcpy( z_inter, z, (MAXN+1)*sizeof(float) );
                                    free_vector( x, 0, MAXN );
                                    free_vector( y, 0, MAXN );
                                    free_vector( z, 0, MAXN );
                                    x = x_inter;
                                    y = y_inter;
                                    z = z_inter;
                                    MAXN += STEP;
				}
		}

	}

        if ( N > 20000 && DRAW_HISTOGRAM == NO )
        {   
            dots = 1;
        }

	if ( LOGLOG == YES )
	  {
	    for ( i=0 ; i < N ; i++ )
	      {
	        if ( y[i] <= 0 || x[i] <= 0 || z[i] <= 0 )
			{
				printf("\033[31m\033[1mLogarithm requested but negative values present ? Abort.\033[0m\n");
				myexit(1);
			}
	        y[i] = log( y[i] );
	        x[i] = log( x[i] );
	        z[i] = log( z[i] );
	      }
          }
	if ( LOG == YES )
	  {
	    for ( i=0 ; i < N ; i++ )
	      {
	        if ( y[i] <= 0 || z[i] <= 0 )
			{
				printf("\033[31m\033[1mLogarithm requested but negative values present ? Abort.\033[0m\n");
				myexit(1);
			}
	        y[i] = log( y[i] );
	        z[i] = log( z[i] );
	      }
          }


        COLOR_DOTS = YES;
        for ( i=0 ; i < N ; i++ )
          {
            if ( z[i] != (int)(z[i]) )
              {
                COLOR_DOTS = NO;
                break;
              }
          }



        if ( COLOR_DOTS == YES )
          {
            DRAW_DENSITY = YES;
          }
        else  
	    three_columns();
	    
	}



	if ( DRAW_DENSITY == YES )
	{
		draw_density();
	}	



	/* XY plots start and end here */
	if ( DRAW_DENSITY == NO && ( columns == 1 || columns == 2 || columns == -2 ))
	{
	
	while ( (dev = qread(&data) ) ) 
		{
			if ( dev == REDRAW )
				{
					if ( FIRST_REDRAW == NO )
					{
					reshapeviewport();
					do_plot_xy();
					}
					else
						FIRST_REDRAW = NO;
				}
			if ( dev == QKEY )
				myexit( 0 );
			if ( dev == PKEY )
			      {
				gl2ppm("plot.ppm");
				if ( HAVE_DONE_KZ == YES && DRAW_KZ == YES )
				  {
				    FILE *kzout;
				    kzout = fopen("Kolmogorov_Zurbenko.dat", "w" );
				    if ( kzout != NULL )
				    {
                                      fprintf(kzout,"# Kolmogorov-Zurbenko filtering with a half-width of %d\n", halfwidth );
				      for ( i=0 ; i < N ; i++ )
				        fprintf(kzout,"%+15.6f %15.6f\n", x[i], kz[i] );
                                      fclose(kzout);
				    }
				  }

				if ( HAVE_DFT != 0 )
				  {
				    FILE *dftout;
				    dftout = fopen("DFT.dat", "w" );
				    if ( dftout != NULL )
				    {
                                      fprintf(dftout,"# Magnitudes of discrete Fourier transform\n" );
				      for ( i=0 ; i < N ; i++ )
				        {
				        if ( HAVE_DFT == 1 )
				           fprintf(dftout,"%16.12f\n", data_DFT1[i] );
                                        else
				           fprintf(dftout,"%16.12f\n", data_DFT2[i] );
				        }
                                      fclose(dftout);
				    }
				  }
                              }
				  
			if ( dev == DKEY )
				{
					dots = dots ^ 1;
					reshapeviewport();
					do_plot_xy();
					usleep(500000);
					qreset();
				}
			if ( dev == TKEY )
				{
				        if ( HAVE_DFT == 0 )
				            {
				              printf("Calculating discete Fourier transform ...\n");
				              start_timer = time(NULL);
				              DFT( N, y);
				              end_timer = time(NULL);
				              if ( (end_timer - start_timer) > 5 )
  				                  printf("Done in %d seconds.\n", (int)(end_timer - start_timer));
                                              else
				                  printf("Done.\n");
                                            }
                                        else
                                              DFT( N, y);
                                              
                                        two_columns();
					reshapeviewport();
					do_plot_xy();
					usleep(500000);
					qreset();
				}
                        if ( dev == PAGEUPKEY )
                                {
					Screencoord	left, right, bottom, top;
                                        float		pscale;
                                        
					getviewport( &left, &right, &bottom, &top );
					if ( fabs(right-left) / 5.0 > fabs(top-bottom) / 5.0 )
					  pscale = fabs(right-left) / 5.0;
                                        else
                                          pscale = fabs(top-bottom) / 5.0;
                                          
                                        minx = minx + (maxx-minx)/pscale;
                                        maxx = maxx - (maxx-minx)/pscale;
                                        miny = miny + (maxy-miny)/pscale;
                                        maxy = maxy - (maxy-miny)/pscale;

					dx = maxx - minx;
					dy = maxy - miny;

					ortho2( minx, maxx, miny, maxy);
					reshapeviewport();
					do_plot_xy();
					qreset();                                  
                                }
                        if ( dev == PAGEDOWNKEY )
                                {
					Screencoord	left, right, bottom, top;
                                        float		pscale;
                                        
					getviewport( &left, &right, &bottom, &top );
					if ( fabs(right-left) / 5.0 > fabs(top-bottom) / 5.0 )
					  pscale = fabs(right-left) / 5.0;
                                        else
                                          pscale = fabs(top-bottom) / 5.0;
                                        minx = minx - (maxx-minx)/pscale;
                                        maxx = maxx + (maxx-minx)/pscale;
                                        miny = miny - (maxy-miny)/pscale;
                                        maxy = maxy + (maxy-miny)/pscale;

					dx = maxx - minx;
					dy = maxy - miny;

					ortho2( minx, maxx, miny, maxy);
					reshapeviewport();
					do_plot_xy();
					qreset();                                  
                                }
                        if ( dev == RIGHTARROWKEY )
                                {
					Screencoord	left, right, bottom, top;
                                        float		pscale;
                                        
					getviewport( &left, &right, &bottom, &top );
					pscale = fabs(right-left) / 5.0;

                                        minx = minx + (maxx-minx)/pscale;
                                        maxx = maxx + (maxx-minx)/pscale;
					ortho2( minx, maxx, miny, maxy);
					do_plot_xy();
					qreset();                                  
                                }
                        if ( dev == LEFTARROWKEY )
                                {
					Screencoord	left, right, bottom, top;
                                        float		pscale;
                                        
					getviewport( &left, &right, &bottom, &top );
					pscale = fabs(right-left) / 5.0;

                                        minx = minx - (maxx-minx)/pscale;
                                        maxx = maxx - (maxx-minx)/pscale;
					ortho2( minx, maxx, miny, maxy);
					do_plot_xy();
					qreset();                                  
                                }
                        if ( dev == DOWNARROWKEY )
                                {
					Screencoord	left, right, bottom, top;
                                        float		pscale;
                                        
					getviewport( &left, &right, &bottom, &top );
					pscale = fabs(top-bottom) / 5.0;

                                        miny = miny - (maxy-miny)/pscale;
                                        maxy = maxy - (maxy-miny)/pscale;
					ortho2( minx, maxx, miny, maxy);
					do_plot_xy();
					qreset();                                  
                                }
                        if ( dev == UPARROWKEY )
                                {
					Screencoord	left, right, bottom, top;
                                        float		pscale;
                                        
					getviewport( &left, &right, &bottom, &top );
					pscale = fabs(top-bottom) / 5.0;

                                
                                        miny = miny + (maxy-miny)/pscale;
                                        maxy = maxy + (maxy-miny)/pscale;
					ortho2( minx, maxx, miny, maxy);
					do_plot_xy();
					qreset();                                  
                                }
			if ( dev == SKEY )
				{
				        if ( HAVE_DONE_KZ == NO )
				          kzfilter();
                                        if ( DRAW_KZ == NO )
                                          DRAW_KZ = YES;
                                        else
					  {
                                          DRAW_KZ = NO;
					  wintitle( title );
                                          }
					reshapeviewport();
					do_plot_xy();
					usleep(500000);
					qreset();
				}
			if ( dev == ZKEY )
				{
					
					HalfWidthOffset += (int)(halfwidth*0.050+0.50)>1?(int)(halfwidth*0.050+0.50):1;
					HAVE_DONE_KZ = NO;
					kzfilter();
					DRAW_KZ = YES;
					reshapeviewport();
					do_plot_xy();
					usleep(5000);
					qreset();
				}
			if ( dev == XKEY )
				{
					HalfWidthOffset -= (int)(halfwidth*0.050+0.50)>1?(int)(halfwidth*0.050+0.50):1;
					HAVE_DONE_KZ = NO;
					kzfilter();
					DRAW_KZ = YES;
					reshapeviewport();
					usleep(5000);
					do_plot_xy();

				}
			if ( dev == FKEY )
				{
					filled_dots = filled_dots ^ 1;
					reshapeviewport();
					do_plot_xy();
					usleep(500000);
					qreset();
				}
			if ( dev == LKEY )
				{
				        PLOT_LABELS = PLOT_LABELS ^ 1;
					reshapeviewport();
					do_plot_xy();
					usleep(500000);
					qreset();
				}
			if ( dev == HKEY )
				{
				        PLOT_HIST = PLOT_HIST ^ 1;
					reshapeviewport();
					do_plot_xy();
					usleep(500000);
					qreset();
				}
			if ( dev == LEFTMOUSE )
				{
					Screencoord	left, right, bottom, top;
					float	xpos, ypos;
					Icoord	orix, oriy;
					
					xpos = getvaluator(MOUSEX);
					ypos = getvaluator(MOUSEY);
					getorigin( &orix, &oriy );
					getviewport( &left, &right, &bottom, &top );
					
					printf("% 15.7f % 15.7f\n", dx*(xpos-orix)/(right-left) + minx, dy*(ypos-oriy)/(top-bottom) + miny );
					usleep(500000);
					qreset();
				}
			if ( dev == RIGHTMOUSE )
				{
					Screencoord	left, right, bottom, top;
					float	xpos1, ypos1;
					float	xpos2, ypos2;
					Icoord	orix, oriy;
					float	x1, x2, y1, y2, inter;
					
					xpos1 = getvaluator(MOUSEX);
					ypos1 = getvaluator(MOUSEY);
					getorigin( &orix, &oriy );
					getviewport( &left, &right, &bottom, &top );
					
					x1 = dx*(xpos1-orix)/(right-left) + minx ;
					y1 = dy*(ypos1-oriy)/(top-bottom) + miny ;
					
					printf("Right-click again to define the other corner ...\n");
					usleep(500000);

					qreset();
					while ( (dev = qread(&data) ) )
						if ( dev == RIGHTMOUSE )
							break;
					
					xpos2 = getvaluator(MOUSEX);
					ypos2 = getvaluator(MOUSEY);
					getorigin( &orix, &oriy );
					getviewport( &left, &right, &bottom, &top );
					
					x2 = dx*(xpos2-orix)/(right-left) + minx ;
					y2 = dy*(ypos2-oriy)/(top-bottom) + miny ;
					
					if ( x1 > x2 )
						{
							inter = x2;
							x2 = x1;
							x1 = inter;
						}
					if ( y1 > y2 )
						{
							inter = y2;
							y2 = y1;
							y1 = inter;
						}	
					
					if ( fabs(xpos2-xpos1) < 2 && fabs(ypos2-ypos1) < 2 )
						{
							two_columns();
						}
					else
					{
					minx = x1;
					maxx = x2;
					miny = y1;
					maxy = y2;                					

					dx = maxx - minx;
					dy = maxy - miny;

					minx = minx - 0.050 * dx;
					maxx = maxx + 0.050 * dx;
					miny = miny - 0.050 * dy;
					maxy = maxy + 0.050 * dy;

					dx = maxx - minx;
					dy = maxy - miny;

					ortho2( minx, maxx, miny, maxy);
					}

					reshapeviewport();
					do_plot_xy();
					usleep(500000);
					qreset();
				}
		}
	}

	/* scatter plots start and end here */
	if ( DRAW_DENSITY == YES )
	{
	while ( (dev = qread(&data) ) ) 
		{
			if ( dev == REDRAW )
				{				
					if ( FIRST_REDRAW == NO )
					{
					reshapeviewport();
					do_draw_density();
					}
					else
						FIRST_REDRAW = NO;
				}
			if ( dev == QKEY )
				myexit( 0 );
			if ( dev == PKEY )
				gl2ppm("plot.ppm");
			if ( dev == DKEY )
				{
					dots = dots ^ 1;
					reshapeviewport();
					do_draw_density();
					usleep(500000);
					qreset();
				}
			if ( dev == LKEY )
				{
				        PLOT_LABELS = PLOT_LABELS ^ 1;
					reshapeviewport();
					do_draw_density();
					usleep(500000);
					qreset();
				}
			if ( dev == FKEY )
				{
					filled_dots = filled_dots ^ 1;
					reshapeviewport();
					do_draw_density();
					usleep(500000);
					qreset();
				}
			if ( dev == EQUALKEY )
				{
				        NOW_PLAYING += STEP_PLAYING;

				        if ( TIMING == NO )
				          {
				            start_timer = time(NULL);
				            TIMING = YES;
				          }
                                        
                                        if ( TIMING == YES )
                                        {
                                          if ( NOW_PLAYING % 100 == 0 && NOW_PLAYING > 0)
                                            {
                                              end_timer = time(NULL);
                                              if ( (int)(end_timer-start_timer) < 5 )
                                                {
                                                  start_timer = end_timer;
                                                  STEP_PLAYING++;
                                                  if ( STEP_PLAYING > 10 )
                                                    {
                                                      STEP_PLAYING = 10;
                                                      TIMING = -1;
                                                    }
                                                }
                                              else
                                                {
                                                  STEP_PLAYING = 1;
                                                  TIMING = NO;
                                                }
                                            }
                                        }
                                        
					reshapeviewport();
					do_draw_density();
					qreset();
				}
			if ( dev == MINUSKEY )
				{
				        NOW_PLAYING--;
				        STEP_PLAYING = 1;
				        TIMING = NO;
				        if ( NOW_PLAYING <= -1 )
				          NOW_PLAYING = N-1;
					reshapeviewport();
					do_draw_density();
					qreset();
				}
                        if ( dev == PAGEUPKEY )
                                {
					Screencoord	left, right, bottom, top;
                                        float		pscale;
                                        
					getviewport( &left, &right, &bottom, &top );
					if ( fabs(right-left) / 5.0 > fabs(top-bottom) / 5.0 )
					  pscale = fabs(right-left) / 5.0;
                                        else
                                          pscale = fabs(top-bottom) / 5.0;
                                        minx = minx + (maxx-minx)/pscale;
                                        maxx = maxx - (maxx-minx)/pscale;
                                        miny = miny + (maxy-miny)/pscale;
                                        maxy = maxy - (maxy-miny)/pscale;

					dx = maxx - minx;
					dy = maxy - miny;

					ortho2( minx, maxx, miny, maxy);
					reshapeviewport();
					do_draw_density();
					qreset();                                  
                                }
                        if ( dev == PAGEDOWNKEY )
                                {
					Screencoord	left, right, bottom, top;
                                        float		pscale;
                                        
					getviewport( &left, &right, &bottom, &top );
					if ( fabs(right-left) / 5.0 > fabs(top-bottom) / 5.0 )
					  pscale = fabs(right-left) / 5.0;
                                        else
                                          pscale = fabs(top-bottom) / 5.0;
                                        minx = minx - (maxx-minx)/pscale;
                                        maxx = maxx + (maxx-minx)/pscale;
                                        miny = miny - (maxy-miny)/pscale;
                                        maxy = maxy + (maxy-miny)/pscale;

					dx = maxx - minx;
					dy = maxy - miny;
					
					ortho2( minx, maxx, miny, maxy);
					reshapeviewport();
					do_draw_density();
					qreset();                                  
                                }
                        if ( dev == RIGHTARROWKEY )
                                {
					Screencoord	left, right, bottom, top;
                                        float		pscale;
                                        
					getviewport( &left, &right, &bottom, &top );
					pscale = fabs(right-left) / 5.0;

                                        minx = minx + (maxx-minx)/pscale;
                                        maxx = maxx + (maxx-minx)/pscale;
					ortho2( minx, maxx, miny, maxy);
					do_draw_density();
					qreset();                                  
                                }
                        if ( dev == LEFTARROWKEY )
                                {
					Screencoord	left, right, bottom, top;
                                        float		pscale;
                                        
					getviewport( &left, &right, &bottom, &top );
					pscale = fabs(right-left) / 5.0;

                                        minx = minx - (maxx-minx)/pscale;
                                        maxx = maxx - (maxx-minx)/pscale;
					ortho2( minx, maxx, miny, maxy);
					do_draw_density();
					qreset();                                  
                                }
                        if ( dev == DOWNARROWKEY )
                                {
					Screencoord	left, right, bottom, top;
                                        float		pscale;
                                        
					getviewport( &left, &right, &bottom, &top );
					pscale = fabs(top-bottom) / 5.0;

                                        miny = miny - (maxy-miny)/pscale;
                                        maxy = maxy - (maxy-miny)/pscale;
					ortho2( minx, maxx, miny, maxy);
					do_draw_density();
					qreset();                                  
                                }
                        if ( dev == UPARROWKEY )
                                {
					Screencoord	left, right, bottom, top;
                                        float		pscale;
                                        
					getviewport( &left, &right, &bottom, &top );
					pscale = fabs(top-bottom) / 5.0;

                                
                                        miny = miny + (maxy-miny)/pscale;
                                        maxy = maxy + (maxy-miny)/pscale;
					ortho2( minx, maxx, miny, maxy);
                                        do_draw_density();
					qreset();                                  
                                }
			if ( dev == LEFTMOUSE )
				{
					Screencoord	left, right, bottom, top;
					float	xpos, ypos;
					Icoord	orix, oriy;
					
					xpos = getvaluator(MOUSEX);
					ypos = getvaluator(MOUSEY);
					getorigin( &orix, &oriy );
					getviewport( &left, &right, &bottom, &top );
					
					printf("% 15.7f % 15.7f\n", dx*(xpos-orix)/(right-left) + minx, dy*(ypos-oriy)/(top-bottom) + miny );
					usleep(500000);
					qreset();
				}
			if ( dev == RIGHTMOUSE )
				{
					Screencoord	left, right, bottom, top;
					float	xpos1, ypos1;
					float	xpos2, ypos2;
					Icoord	orix, oriy;
					float	x1, x2, y1, y2, inter;
					
					xpos1 = getvaluator(MOUSEX);
					ypos1 = getvaluator(MOUSEY);
					getorigin( &orix, &oriy );
					getviewport( &left, &right, &bottom, &top );
					
					x1 = dx*(xpos1-orix)/(right-left) + minx ;
					y1 = dy*(ypos1-oriy)/(top-bottom) + miny ;
					
					printf("Right-click again to define the other corner ...\n");
					usleep(500000);

					qreset();
					while ( (dev = qread(&data) ) )
						if ( dev == RIGHTMOUSE )
							break;
					
					xpos2 = getvaluator(MOUSEX);
					ypos2 = getvaluator(MOUSEY);
					getorigin( &orix, &oriy );
					getviewport( &left, &right, &bottom, &top );
					
					x2 = dx*(xpos2-orix)/(right-left) + minx ;
					y2 = dy*(ypos2-oriy)/(top-bottom) + miny ;
					
					if ( x1 > x2 )
						{
							inter = x2;
							x2 = x1;
							x1 = inter;
						}
					if ( y1 > y2 )
						{
							inter = y2;
							y2 = y1;
							y1 = inter;
						}	
					
					if ( fabs(xpos2-xpos1) < 2 && fabs(ypos2-ypos1) < 2 )
						{
							draw_density();
						}
					else
					{
					minx = x1;
					maxx = x2;
					miny = y1;
					maxy = y2;                					

					dx = maxx - minx;
					dy = maxy - miny;

					minx = minx - 0.050 * dx;
					maxx = maxx + 0.050 * dx;
					miny = miny - 0.050 * dy;
					maxy = maxy + 0.050 * dy;

					dx = maxx - minx;
					dy = maxy - miny;

					ortho2( minx, maxx, miny, maxy);
					}

					reshapeviewport();
					do_draw_density();
					usleep(500000);
					qreset();
				}
			
		}
	}


	/* Three columns are here */
	if ( columns == 3)
	{
	while ( (dev = qread(&data) ) ) 
		{
			if ( dev == REDRAW )
				{
					if ( FIRST_REDRAW == NO )
					{
					reshapeviewport();
					do_plot_two();
					}
					else
						FIRST_REDRAW = NO;
				}
			if ( dev == QKEY )
				myexit( 0 );
			if ( dev == PKEY )
				gl2ppm("plot.ppm");
			if ( dev == DKEY )
				{
					dots = dots ^ 1;
					reshapeviewport();
					do_plot_two();
					usleep(500000);
					qreset();
				}
			if ( dev == FKEY )
				{
					filled_dots = filled_dots ^ 1;
					reshapeviewport();
					do_plot_two();
					usleep(500000);
					qreset();
				}
			if ( dev == LKEY )
				{
				        PLOT_LABELS = PLOT_LABELS ^ 1;
					reshapeviewport();
					do_plot_two();
					usleep(500000);
					qreset();
				}
                        if ( dev == PAGEUPKEY )
                                {
					Screencoord	left, right, bottom, top;
                                        float		pscale;
                                        
					getviewport( &left, &right, &bottom, &top );
					if ( fabs(right-left) / 5.0 > fabs(top-bottom) / 5.0 )
					  pscale = fabs(right-left) / 5.0;
                                        else
                                          pscale = fabs(top-bottom) / 5.0;
                                        minx = minx + (maxx-minx)/pscale;
                                        maxx = maxx - (maxx-minx)/pscale;
                                        miny = miny + (maxy-miny)/pscale;
                                        maxy = maxy - (maxy-miny)/pscale;

					dx = maxx - minx;
					dy = maxy - miny;

					ortho2( minx, maxx, miny, maxy);
					reshapeviewport();
					do_plot_two();
					qreset();                                  
                                }
                        if ( dev == PAGEDOWNKEY )
                                {
					Screencoord	left, right, bottom, top;
                                        float		pscale;
                                        
					getviewport( &left, &right, &bottom, &top );
					if ( fabs(right-left) / 5.0 > fabs(top-bottom) / 5.0 )
					  pscale = fabs(right-left) / 5.0;
                                        else
                                          pscale = fabs(top-bottom) / 5.0;
                                        minx = minx - (maxx-minx)/pscale;
                                        maxx = maxx + (maxx-minx)/pscale;
                                        miny = miny - (maxy-miny)/pscale;
                                        maxy = maxy + (maxy-miny)/pscale;

					dx = maxx - minx;
					dy = maxy - miny;
					
					ortho2( minx, maxx, miny, maxy);
					reshapeviewport();
                                        do_plot_two();
					qreset();                                  
                                }
                        if ( dev == RIGHTARROWKEY )
                                {
					Screencoord	left, right, bottom, top;
                                        float		pscale;
                                        
					getviewport( &left, &right, &bottom, &top );
					pscale = fabs(right-left) / 5.0;

                                        minx = minx + (maxx-minx)/pscale;
                                        maxx = maxx + (maxx-minx)/pscale;
					ortho2( minx, maxx, miny, maxy);
					do_plot_two();
					qreset();                                  
                                }
                        if ( dev == LEFTARROWKEY )
                                {
					Screencoord	left, right, bottom, top;
                                        float		pscale;
                                        
					getviewport( &left, &right, &bottom, &top );
					pscale = fabs(right-left) / 5.0;

                                        minx = minx - (maxx-minx)/pscale;
                                        maxx = maxx - (maxx-minx)/pscale;
					ortho2( minx, maxx, miny, maxy);
					do_plot_two();
					qreset();                                  
                                }
                        if ( dev == DOWNARROWKEY )
                                {
					Screencoord	left, right, bottom, top;
                                        float		pscale;
                                        
					getviewport( &left, &right, &bottom, &top );
					pscale = fabs(top-bottom) / 5.0;

                                        miny = miny - (maxy-miny)/pscale;
                                        maxy = maxy - (maxy-miny)/pscale;
					ortho2( minx, maxx, miny, maxy);
					do_plot_two();
					qreset();                                  
                                }
                        if ( dev == UPARROWKEY )
                                {
					Screencoord	left, right, bottom, top;
                                        float		pscale;
                                        
					getviewport( &left, &right, &bottom, &top );
					pscale = fabs(top-bottom) / 5.0;

                                
                                        miny = miny + (maxy-miny)/pscale;
                                        maxy = maxy + (maxy-miny)/pscale;
					ortho2( minx, maxx, miny, maxy);
					do_plot_two();
					qreset();                                  
                                }
			if ( dev == LEFTMOUSE )
				{
					Screencoord	left, right, bottom, top;
					float	xpos, ypos;
					Icoord	orix, oriy;
					
					xpos = getvaluator(MOUSEX);
					ypos = getvaluator(MOUSEY);
					getorigin( &orix, &oriy );
					getviewport( &left, &right, &bottom, &top );
					
					printf("% 15.7f % 15.7f\n", dx*(xpos-orix)/(right-left) + minx, dy*(ypos-oriy)/(top-bottom) + miny );
					usleep(500000);
					qreset();
				}
			if ( dev == RIGHTMOUSE )
				{
					Screencoord	left, right, bottom, top;
					float	xpos1, ypos1;
					float	xpos2, ypos2;
					Icoord	orix, oriy;
					float	x1, x2, y1, y2, inter;
					
					xpos1 = getvaluator(MOUSEX);
					ypos1 = getvaluator(MOUSEY);
					getorigin( &orix, &oriy );
					getviewport( &left, &right, &bottom, &top );
					
					x1 = dx*(xpos1-orix)/(right-left) + minx ;
					y1 = dy*(ypos1-oriy)/(top-bottom) + miny ;
					
					printf("Right-click again to define the other corner ...\n");
					usleep(500000);

					qreset();
					while ( (dev = qread(&data) ) )
						if ( dev == RIGHTMOUSE )
							break;
					
					xpos2 = getvaluator(MOUSEX);
					ypos2 = getvaluator(MOUSEY);
					getorigin( &orix, &oriy );
					getviewport( &left, &right, &bottom, &top );
					
					x2 = dx*(xpos2-orix)/(right-left) + minx ;
					y2 = dy*(ypos2-oriy)/(top-bottom) + miny ;
					
					if ( x1 > x2 )
						{
							inter = x2;
							x2 = x1;
							x1 = inter;
						}
					if ( y1 > y2 )
						{
							inter = y2;
							y2 = y1;
							y1 = inter;
						}	
					
					if ( fabs(xpos2-xpos1) < 2 && fabs(ypos2-ypos1) < 2 )
						{
							three_columns();
						}
					else
					{
					minx = x1;
					maxx = x2;
					miny = y1;
					maxy = y2;                					

					dx = maxx - minx;
					dy = maxy - miny;

					minx = minx - 0.050 * dx;
					maxx = maxx + 0.050 * dx;
					miny = miny - 0.050 * dy;
					maxy = maxy + 0.050 * dy;

					dx = maxx - minx;
					dy = maxy - miny;

					ortho2( minx, maxx, miny, maxy);
					}

					reshapeviewport();
					do_plot_two();
					usleep(500000);
					qreset();
				}
		}
	}
	
	
	return(1);
}








void	two_columns()
{


        if ( DATA_SET >= 1 )
        {
            int current = 1;
            int starting;
            int finishing;

            i = 1;
            while( current < DATA_SET )
            {
                while( x[i] > x[i-1] && i < N )
                    i++;

                if ( i == N )
                    break;
                i++;
                current++;
            }


            if ( i == N )
            {
                if ( current != DATA_SET )
                    printf("\033[31m\033[1m\nCaution: requested data set not found! Will plot all data.\n\n\033[0m");
            }
            else
            {
                    starting = i-1;
                    while( x[i] > x[i-1] && i < N )
                        i++;
                    finishing = i;
                    for ( i=starting ; i < finishing ; i++ )
                    {
                        x[i-starting] = x[i];
                        y[i-starting] = y[i];
                    }
                    N = finishing - starting;
            }

        }

	for ( i=1 ; i < N ; i++ )
		if ( x[i] < x[i-1] )
			{
				DRAW_DENSITY = YES;
                                dots = YES;
				return;
			}
			

	maxx = x[0];
	minx = maxx;
	maxy = y[0];
	miny = maxy;
	for ( i=1 ; i < N ; i++ )
		{
	
			if ( x[i] > maxx )
				maxx = x[i];
			if ( x[i] < minx )
				minx = x[i];
			if ( y[i] > maxy )
				maxy = y[i];
			if ( y[i] < miny )
				miny = y[i];
		}

	if ( N < 2 )
		{
                        printf("\033[31m\033[1mOnly one point ? Abort.\033[0m\n");
                        myexit( 1 );
		}
	if ( maxx == minx )
		{
			maxx += 0.050 * maxx;
			minx -= 0.050 * minx;
			
			if ( maxx == minx )
			  {
			    maxx =  0.5;
			    minx = -0.5;
			  }

		}
	if ( maxy == miny )
		{
			maxy += 0.050 * maxy;
			miny -= 0.050 * miny;

			if ( maxy == miny )
			  {
			    maxy =  0.5;
			    miny = -0.5;
			  }
		}

	if ( HAVE_MIN_MAX == YES && DRAW_HISTOGRAM == NO )
	  {
	    miny = MIN_FROM_CML;
	    maxy = MAX_FROM_CML;
	  }

	dx = maxx - minx;
	dy = maxy - miny;

	minx = minx - 0.050 * dx;
	maxx = maxx + 0.050 * dx;
	miny = miny - 0.050 * dy;
	maxy = maxy + 0.050 * dy;

	dx = maxx - minx;
	dy = maxy - miny;

	sprintf( title, "Limits are %+6.4E to %+6.4E on x, %+6.4E to %+6.4E on y", minx, maxx, miny, maxy );


	if ( FIRST == 0 )
	{
        minsize( MINX_SIZE, MINY_SIZE );
        winid = winopen( title );
        if ( winid == -1 )
                {
                        printf("\033[37m\033[1mCan't open graphics window. Abort.\033[0m\n");
                        myexit( 1 );
		}
		
        winset ( winid );
        FIRST = 1;
        }

	ortho2( minx, maxx, miny, maxy);
	doublebuffer();
	gconfig();
	color( 0 );
	clear();
        if ( LARGE_LABELS == 0 )
            loadXfont(4711, "fixed");
        else
            loadXfont(4711, "10x20");
        font(4711);
        color(7);
        cmov2( minx, miny );
/*        charstr( "Data loaded. Rendering ..." );
*/        gflush();
        gsync();
	swapbuffers();
	
	qdevice(QKEY);
	qdevice(DKEY);
	qdevice(FKEY);
	qdevice(LKEY);
	qdevice(PKEY);
	qdevice(SKEY);
	qdevice(HKEY);
	qdevice(ZKEY);
	qdevice(XKEY);	
	qdevice(TKEY);	
        qdevice(PAGEUPKEY);	
        qdevice(PAGEDOWNKEY);	
	qdevice(DOWNARROWKEY);	
	qdevice(UPARROWKEY);	
	qdevice(RIGHTARROWKEY);	
	qdevice(LEFTARROWKEY);	
	qdevice(LEFTMOUSE);
	qdevice(RIGHTMOUSE);
	do_plot_xy();



	return;
}




void 	do_plot_xy()
{
	int i;
	Screencoord     left, right, bottom, top;


	if ( have_plot_labels == 0 && PLOT_LABELS == 1 )
	{
	float	stepx, stepy;
	float 	roundminy, roundminx;
	
        gflush();
        gsync();
        backbuffer( 1 );

	have_plot_labels = 1;
	
	color( 6 );
	mapcolor( 201  ,50, 50, 50 );
	linewidth( 1 );
	deflinestyle( 2 , 0x0101 );
	setlinestyle( 2 );

	stepx = pow( 10.0 , floor( log10( (maxx-minx)/5 )) );
	stepy = pow( 10.0 , floor( log10( (maxy-miny)/5 )) );
	
	stepx *= (int)( ( (maxx-minx) / stepx ) / 5 );
	stepy *= (int)( ( (maxy-miny) / stepy ) / 5 );
	
	roundminx = stepx * (int)(minx / stepx);
	roundminy = stepy * (int)(miny / stepy);
	
	for ( i=0 ; i < (int)((maxy-miny)/stepy+2.5)  ; i++ )
	  {
	    if ( i==0 && roundminy < miny ) 
	    		;
	    else
	    {
	    if ( PLOT_LIGHTGRID == 0 )
	    {
	    move2( minx, i*stepy+roundminy );
	    draw2( maxx, i*stepy+roundminy );
	    }
	    else
	    {
	    color( 201 );
	    setlinestyle( 0 );
	    move2( minx, i*stepy+roundminy );
	    draw2( maxx, i*stepy+roundminy );
	    setlinestyle( 2 );
	    color( 6 );
	    }
	    
	    if ( i*stepy+roundminy < 0.9850*(maxy-miny)+miny )
	    {
	    cmov2( minx, i*stepy+roundminy );
	    if ( fabs( maxy-miny ) > 9999999 || fabs(maxy) > 9999999 || fabs(miny) > 9999999 )
	    	sprintf( label, "%+6.4E", i*stepy+roundminy );
	    else if ( fabs( maxy-miny ) > 50 )
	    	sprintf( label, "%.0f", i*stepy+roundminy );
	    else if ( fabs( maxy-miny ) > 5 )
	    	sprintf( label, "%.1f", i*stepy+roundminy );
	    else if ( fabs( maxy-miny ) > 1 ) 
	        sprintf( label, "%.2f", i*stepy+roundminy );
	    else
	      	sprintf( label, "%+6.4E", i*stepy+roundminy );
	    charstr( label );
	    }

	    }
	  }
	for ( i=0 ; i < (int)((maxx-minx)/stepx+2.5) ; i++ )
	  {
	    if ( i==0 && roundminx < minx )
	    		;
	    else
	    {
	    if ( PLOT_LIGHTGRID == 0 )
	    {
	    move2( i*stepx+roundminx, miny );
	    draw2( i*stepx+roundminx, maxy );
	    }
	    else
	    {
	    color( 201 );
	    setlinestyle( 0 );
	    move2( i*stepx+roundminx, miny );
	    draw2( i*stepx+roundminx, maxy );
	    setlinestyle( 2 );
	    color( 6 );
	    }
	    
	    if ( i*stepx+roundminx < 0.9750*(maxx-minx)+minx )
	    {
	    cmov2( i*stepx+roundminx, miny );
	    if ( fabs( maxx-minx ) > 9999999 || fabs(maxx) > 9999999 || fabs(minx) > 9999999 )
	    	sprintf( label, "%+6.4E", i*stepx+roundminx );
	    else if ( fabs( maxx-minx ) > 50 )
	    	sprintf( label, "%.0f", i*stepx+roundminx );
	    else if ( fabs( maxx-minx ) > 5 )
	    	sprintf( label, "%.1f", i*stepx+roundminx );
	    else if ( fabs( maxx-minx ) > 1 ) 
	        sprintf( label, "%.2f", i*stepx+roundminx );
	    else
	      	sprintf( label, "%+6.4E", i*stepx+roundminx );
	    charstr( label );
	    }
	    }
	  }
        setlinestyle( 0 );

	color( 7 );
	linewidth( 1 );
	if ( minx < 0 && maxx > 0 )
		{
			move2( 0.0, miny );
			draw2( 0.0, maxy );
		}

	if ( miny < 0 && maxy > 0 )
		{
			move2( minx, 0.0 );
			draw2( maxx, 0.0 );
		}
		

        gflush();
        gsync();
        swapbuffers();
	return;                			

	}
	


        frontbuffer( 1 );
        if ( COLOR == YES )
          color(0);
        else
          color(7);
        cmov2( minx, miny );
        charstr( "Rendering ..." );
        gflush();
        gsync();
        backbuffer( 1 );
	
        color( 0 );
        clear();
	linewidth( 1 );


	if ( PLOT_LABELS == 1 )
	{
	float	stepx, stepy;
	float 	roundminy, roundminx;
	
	have_plot_labels = 1;
	
	color( 6 );
	mapcolor( 201  ,50, 50, 50 );
	linewidth( 1 );
	deflinestyle( 2 , 0x0101 );
	setlinestyle( 2 );

	stepx = pow( 10.0 , floor( log10( (maxx-minx)/5 )) );
	stepy = pow( 10.0 , floor( log10( (maxy-miny)/5 )) );
	
	stepx *= (int)( ( (maxx-minx) / stepx ) / 5 );
	stepy *= (int)( ( (maxy-miny) / stepy ) / 5 );
	
	roundminx = stepx * (int)(minx / stepx);
	roundminy = stepy * (int)(miny / stepy);
	
	for ( i=0 ; i < (int)((maxy-miny)/stepy+2.5)  ; i++ )
	  {
	    if ( i==0 && roundminy < miny ) 
	    		;
	    else
	    {
	    if ( PLOT_LIGHTGRID == 0 )
	    {	    
	    move2( minx, i*stepy+roundminy );
	    draw2( maxx, i*stepy+roundminy );
	    }
	    else
	    {
	    color( 201 );
	    setlinestyle( 0 );
	    move2( minx, i*stepy+roundminy );
	    draw2( maxx, i*stepy+roundminy );
	    setlinestyle( 2 );
	    color( 6 );
	    }
	    
	    if ( i*stepy+roundminy < 0.9850*(maxy-miny)+miny )
	    {
	    cmov2( minx, i*stepy+roundminy );
	    if ( fabs( maxy-miny ) > 9999999 || fabs(maxy) > 9999999 || fabs(miny) > 9999999 )
	    	sprintf( label, "%+6.4E", i*stepy+roundminy );
	    else if ( fabs( maxy-miny ) > 50 )
	    	sprintf( label, "%.0f", i*stepy+roundminy );
	    else if ( fabs( maxy-miny ) > 5 )
	    	sprintf( label, "%.1f", i*stepy+roundminy );
	    else if ( fabs( maxy-miny ) > 1 ) 
	        sprintf( label, "%.2f", i*stepy+roundminy );
	    else
	      	sprintf( label, "%+6.4E", i*stepy+roundminy );
	    charstr( label );
	    }

	    }
	  }
	for ( i=0 ; i < (int)((maxx-minx)/stepx+2.5) ; i++ )
	  {
	    if ( i==0 && roundminx < minx )
	    		;
	    else
	    {
	    if ( PLOT_LIGHTGRID == 0 )
	    {
	    move2( i*stepx+roundminx, miny );
	    draw2( i*stepx+roundminx, maxy );
	    }
	    else
	    {
	    color( 201 );
	    setlinestyle( 0 );
	    move2( i*stepx+roundminx, miny );
	    draw2( i*stepx+roundminx, maxy );
	    setlinestyle( 2 );
	    color( 6 );
	    }
	    
	    if ( i*stepx+roundminx < 0.9750*(maxx-minx)+minx )
	    {
	    cmov2( i*stepx+roundminx, miny );
	    if ( fabs( maxx-minx ) > 9999999 || fabs(maxx) > 9999999 || fabs(minx) > 9999999 )
	    	sprintf( label, "%+6.4E", i*stepx+roundminx );
	    else if ( fabs( maxx-minx ) > 50 )
	    	sprintf( label, "%.0f", i*stepx+roundminx );
	    else if ( fabs( maxx-minx ) > 5 )
	    	sprintf( label, "%.1f", i*stepx+roundminx );
	    else if ( fabs( maxx-minx ) > 1 ) 
	        sprintf( label, "%.2f", i*stepx+roundminx );
	    else
	      	sprintf( label, "%+6.4E", i*stepx+roundminx );
	    charstr( label );
	    }
	    }
	  }
        setlinestyle( 0 );
        }
        else
		have_plot_labels = 0;
		

	if ( DRAW_HISTOGRAM == YES && PLOT_HIST == YES )
	{
		mapcolor( 201  ,150, 150, 150 );
		color( 201 );
		for ( i=0 ; i < N ; i++ )
		{
			move2( x[i]-hist_step/2, 0 );
			draw2( x[i]-hist_step/2, y[i] );
			draw2( x[i]+hist_step/2, y[i] );
			draw2( x[i]+hist_step/2, 0 );
		}
	}


        mapcolor( 200  ,243, 146, 46 );
	color( 200 );

	if ( dots == NO && filled_dots == NO )
	{
            i = 0;
            while ( x[i] < minx && i < N )
                i++;
       
            if ( i > 0 )
                i--;

	    move2( x[i], y[i] );
	    while ( x[i] < maxx && i < N )
                {
		    draw2( x[i], y[i] );
                    i++;
                }

            if ( i < N )
                draw2( x[i], y[i] );
	}
	else if ( filled_dots == YES )
	{
	    getviewport( &left, &right, &bottom, &top );

            i = 0;
            while ( x[i] < minx && i < N )
                i++;
            
            if ( i > 0 )
                i--;

	    while ( x[i] < maxx && i < N )
                {
		    arcxf( x[i], y[i], 4.0*dx/(right-left), 4.0*dy/(top-bottom), 0, 3600 );
                    i++;
                }

            if ( i < N )
                arcxf( x[i], y[i], 4.0*dx/(right-left), 4.0*dy/(top-bottom), 0, 3600 );
            
            if ( dots == NO ) 
            {
                i = 0;
                while ( x[i] < minx && i < N )
                    i++;
        
                if ( i > 0 )
                    i--;

                move2( x[i], y[i] );
                while ( x[i] < maxx && i < N )
                {
                    draw2( x[i], y[i] );
                    i++;
                }

                if ( i < N )
                    draw2( x[i], y[i] );
            }
	}
	else
	{

            i = 0;
            while ( x[i] < minx && i < N )
                i++;
        
            if ( i > 0 )
                i--;

	    while ( x[i] < maxx && i < N )
                {
		    pnt2( x[i], y[i] );
                    i++;
                }

            if ( i < N )
                pnt2( x[i], y[i] );
	}


	if ( DRAW_KZ == YES )
	  {
	    getviewport( &left, &right, &bottom, &top );
	    color(2);

            
            
            i = 0;
            while ( x[i] < minx && i < N )
                i++;
            
            if ( i > 0 )
                i--;

	    while ( x[i] < maxx && i < N )
                {
                    arcxf( x[i], kz[i], 4.0*dx/(right-left), 4.0*dy/(top-bottom), 0, 3600 );
                    i++;
                }

            if ( i < N )
                arcxf( x[i], kz[i], 4.0*dx/(right-left), 4.0*dy/(top-bottom), 0, 3600 );
            

            if ( N < 500 )
            {
            i = 0;
            while ( x[i] < minx && i < N )
                i++;
       
            if ( i > 0 )
                i--;

	    move2( x[i], kz[i] );
	    while ( x[i] < maxx && i < N )
                {
		    draw2( x[i], kz[i] );
                    i++;
                }

            if ( i < N )
                draw2( x[i], kz[i] );
            }

	  }
	



	color( 7 );
	linewidth( 1 );
	if ( minx < 0 && maxx > 0 )
		{
			move2( 0.0, miny );
			draw2( 0.0, maxy );
		}

	if ( miny < 0 && maxy > 0 )
		{
			move2( minx, 0.0 );
			draw2( maxx, 0.0 );
		}
		


        gflush();
        gsync();
	swapbuffers();

	return;                			
}






void	draw_density()
{
	static int SAVED_MATRIX = NO;
	

	maxx = x[0];
	minx = maxx;
	maxy = y[0];
	miny = maxy;
	for ( i=1 ; i < N ; i++ )
		{
	
			if ( x[i] > maxx )
				maxx = x[i];
			if ( x[i] < minx )
				minx = x[i];
			if ( y[i] > maxy )
				maxy = y[i];
			if ( y[i] < miny )
				miny = y[i];
		}

	if ( N < 2 )
		{
                        printf("\033[31m\033[1mOnly one point ? Abort.\033[0m\n");
                        myexit( 1 );
		}



	if ( maxx == minx )
		{
			maxx += 0.050 * maxx;
			minx -= 0.050 * minx;
		}
	if ( maxy == miny )
		{
			maxy += 0.050 * maxy;
			miny -= 0.050 * miny;
		}
		

	dx = maxx - minx;
	dy = maxy - miny;

	minx = minx - 0.050 * dx;
	maxx = maxx + 0.050 * dx;
	miny = miny - 0.050 * dy;
	maxy = maxy + 0.050 * dy;

	dx = maxx - minx;
	dy = maxy - miny;

	sprintf( title, "Limits are %+6.4E to %+6.4E on x, %+6.4E to %+6.4E on y", minx, maxx, miny, maxy );


	if ( FIRST == 0 )
	{
        minsize( MINS_SIZE, MINS_SIZE );
        keepaspect( MINS_SIZE, MINS_SIZE );
        winid = winopen( title );
        if ( winid == -1 )
                {
                        printf("\033[37m\033[1mCan't open graphics window. Abort.\033[0m\n");
                        myexit( 1 );
		}
		
        winset ( winid );
        FIRST = 1;
        }
        
	ortho2( minx, maxx, miny, maxy);
	doublebuffer();
	gconfig();
	color( 0 );
	clear();

        if ( N < 2000 )
            filled_dots = YES;

        if ( LARGE_LABELS == 0 )
            loadXfont(4711, "fixed");
        else
            loadXfont(4711, "10x20");
        font(4711);
        color(7);
        cmov2( minx, miny );
/*        charstr( "Data loaded. Rendering ..." );
*/        gflush();
        gsync();
	swapbuffers();
	
	qdevice(QKEY);
	qdevice(DKEY);
	qdevice(LKEY);
	qdevice(FKEY);
	qdevice(PKEY);

	qdevice(EQUALKEY);
	qdevice(MINUSKEY);

        qdevice(PAGEUPKEY);	
        qdevice(PAGEDOWNKEY);	
	qdevice(DOWNARROWKEY);	
	qdevice(UPARROWKEY);	
	qdevice(RIGHTARROWKEY);	
	qdevice(LEFTARROWKEY);	

	qdevice(LEFTMOUSE);
	qdevice(RIGHTMOUSE);
	do_draw_density();


	if ( SAVED_MATRIX == NO )
	{
		float	*inter;
		float	**denmat;
		float	max, min;
		int	h_N;
		FILE	*densfile;
		int	k, l;
		float	gmeanx, gmeany;	
		double 	sx=0, sy=0, sx2=0, sy2=0, sxy=0, r;
		float	xx, yy;

		for ( i=0 ; i < N ; i++ )
		{
		xx = x[i];
                yy = y[i];
		sx = sx + xx;
		sy = sy + yy;
		sx2 = sx2 + xx * xx;
		sy2 = sy2 + yy * yy;
		sxy = sxy + xx * yy;
		}

		r = ( N*sxy - sx*sy ) / sqrt( (N*sx2 - sx*sx)*( N*sy2-sy*sy ) );
		printf("Correlation coefficient for %d pairs is %+10.8f\n", N, r);


		
		inter = vector( 0, N-1 );

		gmeanx = 0.0;
		gmeany = 0.0;
		for ( i=0 ; i < N ; i++ )
		{
			gmeanx += x[i];
			gmeany += y[i];
		}
		gmeanx /= N;
		gmeany /= N;
		
                for ( i=0 ; i < N ; i++ )
                       inter[i] = (x[i]-gmeanx)*(x[i]-gmeanx) + (y[i]-gmeany)*(y[i]-gmeany) ;

		sort( N, inter-1 );

		min = inter[0];
		max = inter[N-1];

		if ( FINE_DENSITY_GRID == NO )
			h_N = 2*(int)( sqrt( (max-min) / (2*(inter[ (int)(3.0*N/4.0 + 0.5)]-inter[ (int)(N/4.0 +0.5) ]) / pow( N, 1.0l/3.0l))) + 0.50 );
		else
			h_N = (int)( (max-min) / (2*(inter[ (int)(3.0*N/4.0 + 0.5)]-inter[ (int)(N/4.0 +0.5) ]) / pow( N, 1.0l/3.0l)) + 0.50 );
		

                if ( h_N > 4 * (int)((int)( pow( (double)( N ), 1.0/3.0 ) +0.50) / 2 +0.50) )
			h_N = 4 * (int)((int)( pow( (double)( N ), 1.0/3.0 ) +0.50) / 2 +0.50);

		
		free_vector( inter, 0, N );

		if ( h_N > 3 )
		{
		
		denmat = matrix( 0, h_N, 0, h_N );

		for ( k=0 ; k <= h_N ; k++ )
		for ( l=0 ; l <= h_N ; l++ )
			denmat[k][l] = 0;

		for ( i=0 ; i < N ; i++ )
			denmat[ (int)( h_N*((x[i]-minx)/(maxx-minx)) + 0.5 ) ][ (int)( h_N*((y[i]-miny)/(maxy-miny)) + 0.5 ) ]++;


			
		densfile = fopen( "density.matrix", "w" );
		if ( densfile != NULL )
		{
		for ( k=0 ; k <= h_N ; k++ )
		{
			for ( l=0 ; l <= h_N ; l++ )
				fprintf( densfile, "%8d ", (int)(denmat[l][h_N-k] + 0.50) );
			
			fprintf(densfile, "\n" );
		}
		fclose( densfile );
		}

		densfile = fopen( "density.log.matrix", "w" );
		if ( densfile != NULL )
		{
		for ( k=0 ; k <= h_N ; k++ )
		{
			for ( l=0 ; l <= h_N ; l++ )
				if ( denmat[l][h_N-k] >= 1 )
					fprintf( densfile, "%11.7f ", log(denmat[l][h_N-k]) );
				else
					fprintf( densfile, "%11.7f ", 0.0 );
			
			fprintf(densfile, "\n" );
		}
		fclose( densfile );
		}

		free_matrix( denmat, 0, h_N, 0, h_N );
		
		if ( densfile != NULL )
			fprintf( stderr, "\033[32m\033[1mDensity matrices saved. Plot with '-cc' flag.\033[0m\n" );		
		
		}
		SAVED_MATRIX = YES;

	}



	return;
}




void 	do_draw_density()
{
	int i;
	Screencoord     left, right, bottom, top;


	if ( have_plot_labels == 0 && PLOT_LABELS == 1 )
	{
	float	stepx, stepy;
	float 	roundminy, roundminx;
	
        gflush();
        gsync();
        backbuffer( 1 );

	have_plot_labels = 1;
	
	color( 6 );
	linewidth( 1 );
	mapcolor( 201  ,50, 50, 50 );
	deflinestyle( 2 , 0x0101 );
	setlinestyle( 2 );

	stepx = pow( 10.0 , floor( log10( (maxx-minx)/5 )) );
	stepy = pow( 10.0 , floor( log10( (maxy-miny)/5 )) );
	
	stepx *= (int)( ( (maxx-minx) / stepx ) / 5 );
	stepy *= (int)( ( (maxy-miny) / stepy ) / 5 );
	
	roundminx = stepx * (int)(minx / stepx);
	roundminy = stepy * (int)(miny / stepy);
	
	for ( i=0 ; i < (int)((maxy-miny)/stepy+2.5)  ; i++ )
	  {
	    if ( i==0 && roundminy < miny ) 
	    		;
	    else
	    {
	    if ( PLOT_LIGHTGRID == 0 )
	    {          
	    move2( minx, i*stepy+roundminy );
	    draw2( maxx, i*stepy+roundminy );
	    }
	    else
	    {
	    color( 201 );
	    setlinestyle( 0 );
	    move2( minx, i*stepy+roundminy );
	    draw2( maxx, i*stepy+roundminy );
	    setlinestyle( 2 );
	    color( 6 );
	    }
	    
	    if ( i*stepy+roundminy < 0.9850*(maxy-miny)+miny )
	    {
	    cmov2( minx, i*stepy+roundminy );
	    if ( fabs( maxy-miny ) > 9999999 || fabs(maxy) > 9999999 || fabs(miny) > 9999999 )
	    	sprintf( label, "%+6.4E", i*stepy+roundminy );
	    else if ( fabs( maxy-miny ) > 50 )
	    	sprintf( label, "%.0f", i*stepy+roundminy );
	    else if ( fabs( maxy-miny ) > 5 )
	    	sprintf( label, "%.1f", i*stepy+roundminy );
	    else if ( fabs( maxy-miny ) > 1 ) 
	        sprintf( label, "%.2f", i*stepy+roundminy );
	    else
	      	sprintf( label, "%+6.4E", i*stepy+roundminy );
	    charstr( label );
	    }

	    }

	  }

	for ( i=0 ; i < (int)((maxx-minx)/stepx+2.5) ; i++ )
	  {
	    if ( i==0 && roundminx < minx )
	    		;
	    else
	    {
	    if ( PLOT_LIGHTGRID == 0 )
	    {
	    move2( i*stepx+roundminx, miny );
	    draw2( i*stepx+roundminx, maxy );
	    }
	    else
	    {
	    color( 201 );
	    setlinestyle( 0 );
	    move2( i*stepx+roundminx, miny );
	    draw2( i*stepx+roundminx, maxy );
	    setlinestyle( 2 );
	    color( 6 );
	    }
	    
	    if ( i*stepx+roundminx < 0.9750*(maxx-minx)+minx )
	    {
	    cmov2( i*stepx+roundminx, miny );
	    if ( fabs( maxx-minx ) > 9999999 || fabs(maxx) > 9999999 || fabs(minx) > 9999999 )
	    	sprintf( label, "%+6.4E", i*stepx+roundminx );
	    else if ( fabs( maxx-minx ) > 50 )
	    	sprintf( label, "%.0f", i*stepx+roundminx );
	    else if ( fabs( maxx-minx ) > 5 )
	    	sprintf( label, "%.1f", i*stepx+roundminx );
	    else if ( fabs( maxx-minx ) > 1 ) 
	        sprintf( label, "%.2f", i*stepx+roundminx );
	    else
	      	sprintf( label, "%+6.4E", i*stepx+roundminx );
	    charstr( label );
	    }
	    }
	  }
        setlinestyle( 0 );


	color( 7 );
	linewidth( 2 );
	move2( minx, miny );
	draw2( maxx, miny );
	draw2( maxx, maxy );
	draw2( minx, maxy );
	draw2( minx, miny );


        gflush();
        gsync();
        swapbuffers();
	return;                			

	}


        if ( NOW_PLAYING >= 0 )
          {
            if ( NOW_PLAYING >= N )
              NOW_PLAYING = 0;
          }

        frontbuffer( 1 );
        if ( COLOR == YES )
          color(0);
        else
          color(7);
        cmov2( minx, miny );
        
        charstr( "Rendering ..." );

        gflush();
        gsync();
        backbuffer( 1 );


	
        color( 0 );
        clear();


	if ( PLOT_LABELS == 1 )
	{
	float	stepx, stepy;
	float 	roundminy, roundminx;
	
	have_plot_labels = 1;

	color( 6 );
	linewidth( 1 );
	mapcolor( 201  ,50, 50, 50 );
	deflinestyle( 2 , 0x0101 );
	setlinestyle( 2 );

	stepx = pow( 10.0 , floor( log10( (maxx-minx)/5 )) );
	stepy = pow( 10.0 , floor( log10( (maxy-miny)/5 )) );
	
	stepx *= (int)( ( (maxx-minx) / stepx ) / 5 );
	stepy *= (int)( ( (maxy-miny) / stepy ) / 5 );
	
	roundminx = stepx * (int)(minx / stepx);
	roundminy = stepy * (int)(miny / stepy);
	
	for ( i=0 ; i < (int)((maxy-miny)/stepy+2.5)  ; i++ )
	  {
	    if ( i==0 && roundminy < miny ) 
	    		;
	    else
	    {
	    if ( PLOT_LIGHTGRID == 0 )
	    {	    
	    move2( minx, i*stepy+roundminy );
	    draw2( maxx, i*stepy+roundminy );
            }
            else
            {
            color( 201 );
            setlinestyle( 0 );
	    move2( minx, i*stepy+roundminy );
	    draw2( maxx, i*stepy+roundminy );
	    setlinestyle( 2 );
	    color( 6 );
            }
            
	    if ( i*stepy+roundminy < 0.9850*(maxy-miny)+miny )
	    {
	    cmov2( minx, i*stepy+roundminy );
	    if ( fabs( maxy-miny ) > 9999999 || fabs(maxy) > 9999999 || fabs(miny) > 9999999 )
	    	sprintf( label, "%+6.4E", i*stepy+roundminy );
	    else if ( fabs( maxy-miny ) > 50 )
	    	sprintf( label, "%.0f", i*stepy+roundminy );
	    else if ( fabs( maxy-miny ) > 5 )
	    	sprintf( label, "%.1f", i*stepy+roundminy );
	    else if ( fabs( maxy-miny ) > 1 ) 
	        sprintf( label, "%.2f", i*stepy+roundminy );
	    else
	      	sprintf( label, "%+6.4E", i*stepy+roundminy );
	    charstr( label );
	    }

	    }
	  }
	for ( i=0 ; i < (int)((maxx-minx)/stepx+2.5) ; i++ )
	  {
	    if ( i==0 && roundminx < minx )
	    		;
	    else
	    {
	    if ( PLOT_LIGHTGRID == 0 )
	    {
	    move2( i*stepx+roundminx, miny );
	    draw2( i*stepx+roundminx, maxy );
	    }
	    else
	    {
	    color( 201 );
	    setlinestyle( 0 );
	    move2( i*stepx+roundminx, miny );
	    draw2( i*stepx+roundminx, maxy );
	    setlinestyle( 2 );
	    color( 6 );
	    }
	    
	    if ( i*stepx+roundminx < 0.9750*(maxx-minx)+minx )
	    {
	    cmov2( i*stepx+roundminx, miny );
	    if ( fabs( maxx-minx ) > 9999999 || fabs(maxx) > 9999999 || fabs(minx) > 9999999 )
	    	sprintf( label, "%+6.4E", i*stepx+roundminx );
	    else if ( fabs( maxx-minx ) > 50 )
	    	sprintf( label, "%.0f", i*stepx+roundminx );
	    else if ( fabs( maxx-minx ) > 5 )
	    	sprintf( label, "%.1f", i*stepx+roundminx );
	    else if ( fabs( maxx-minx ) > 1 ) 
	        sprintf( label, "%.2f", i*stepx+roundminx );
	    else
	      	sprintf( label, "%+6.4E", i*stepx+roundminx );
	    charstr( label );
	    }
	    }
	  }
        setlinestyle( 0 );
        }
		else
			have_plot_labels = 0;


        mapcolor( 200  ,243, 146, 46 );
	color( 200 );

        getviewport( &left, &right, &bottom, &top );
        
        if ( dots == NO )
        {

        if ( filled_dots == YES )
            for ( i=0 ; i < N ; i++ )
                arcxf( x[i], y[i], 6.0*dx/(right-left), 6.0*dy/(top-bottom), 0, 3600 );

        linewidth( 1 );
	move2( x[0], y[0] );
	for ( i=1 ; i < N ; i++ )
		draw2( x[i], y[i] );
        }                
	else
	{
	for ( i=0 ; i < N ; i++ )
		{
		  if ( COLOR_DOTS == YES )
		    color( (int)(z[i])%7 + 1 );
		  
                  if ( filled_dots == YES &&  x[i] >= minx && x[i] <= maxx && y[i] >= miny && y[i] <= maxy )
                    arcxf( x[i], y[i], 6.0*dx/(right-left), 6.0*dy/(top-bottom), 0, 3600 );

                  if ( filled_dots == NO && x[i] >= minx && x[i] <= maxx && y[i] >= miny && y[i] <= maxy )
		    pnt2( x[i], y[i] );
                }
	}


	if ( NOW_PLAYING >= 0 )
	{
	getviewport( &left, &right, &bottom, &top );
	color( 6 );
	for ( i=NOW_PLAYING-STEP_PLAYING+1 ; i <= NOW_PLAYING ; i++ )
          arcxf( x[i], y[i], 6.0*dx/(right-left), 6.0*dy/(top-bottom), 0, 3600 );
	}





	color( 7 );
	linewidth( 2 );
	move2( minx, miny );
	draw2( maxx, miny );
	draw2( maxx, maxy );
	draw2( minx, maxy );
	draw2( minx, miny );


        if ( NOW_PLAYING >= 0 )
          {
              
	    sprintf( label, "%10d", NOW_PLAYING+1 );
	    charstr( label );
          }

	
        gflush();
        gsync();
	swapbuffers();

	return;                			
}





void	three_columns()
{


        if ( DATA_SET >= 1 )
        {
            int current = 1;
            int starting;
            int finishing;

            i = 1;
            while( current < DATA_SET )
            {
                while( x[i] > x[i-1] && i < N )
                    i++;

                if ( i == N )
                    break;
                i++;
                current++;
            }


            if ( i == N )
            {
                if ( current != DATA_SET )
                    printf("\033[31m\033[1m\nCaution: requested data set not found! Will plot all data.\n\n\033[0m");
            }
            else
            {
                    starting = i-1;
                    while( x[i] > x[i-1] && i < N )
                        i++;
                    finishing = i;
                    for ( i=starting ; i < finishing ; i++ )
                    {
                        x[i-starting] = x[i];
                        y[i-starting] = y[i];
                        z[i-starting] = z[i];
                    }
                    N = finishing - starting;
            }

        }



	if ( ERROR_BARS == YES )
	{
	maxx = x[0];
	minx = maxx;
	maxy = y[0];
	miny = maxy;
	for ( i=1 ; i < N ; i++ )
		{
	
			if ( x[i] > maxx )
				maxx = x[i];
			if ( x[i] < minx )
				minx = x[i];
			if ( y[i]+z[i] > maxy )
				maxy = y[i]+z[i];
			if ( y[i]-z[i] < miny )
				miny = y[i]-z[i];
		}
	
	}
	else
	{
	if ( AUTOSCALE == NO )
	{
	maxx = x[0];
	minx = maxx;
	maxy = y[0];
	miny = maxy;
	for ( i=1 ; i < N ; i++ )
		{
	
			if ( x[i] > maxx )
				maxx = x[i];
			if ( x[i] < minx )
				minx = x[i];
			if ( y[i] > maxy )
				maxy = y[i];
			if ( y[i] < miny )
				miny = y[i];
			if ( z[i] > maxy )
				maxy = z[i];
			if ( z[i] < miny )
				miny = z[i];
		}
	}
	else
	{				/* Scale/translate second graph to overlap range of first */
	
	maxx = x[0];
	minx = maxx;
	maxy = y[0];
	miny = maxy;
	maxz = z[0];
	minz = maxz;
	for ( i=1 ; i < N ; i++ )
		{
	
			if ( x[i] > maxx )
				maxx = x[i];
			if ( x[i] < minx )
				minx = x[i];
			if ( y[i] > maxy )
				maxy = y[i];
			if ( y[i] < miny )
				miny = y[i];
			if ( z[i] > maxz )
				maxz = z[i];
			if ( z[i] < minz )
				minz = z[i];
		}
	
	scalesecond = (maxy-miny) / (maxz-minz);
	transsecond = miny - minz * scalesecond;
	
	for ( i=0 ; i < N ; i++ )
		z[i] = scalesecond * z[i] + transsecond;
	
	
	}
	
	}



	if ( N < 2 )
		{
                        printf("\033[31m\033[1mOnly one point ? Abort.\033[0m\n");
                        myexit( 1 );
		}
	if ( maxx == minx )
		{
			maxx += 0.050 * maxx;
			minx -= 0.050 * minx;
		}
	if ( maxy == miny )
		{
			maxy += 0.050 * maxy;
			miny -= 0.050 * miny;
		}

	if ( HAVE_MIN_MAX == YES && DRAW_HISTOGRAM == NO )
	  {
	    miny = MIN_FROM_CML;
	    maxy = MAX_FROM_CML;
	  }

	dx = maxx - minx;
	dy = maxy - miny;

	minx = minx - 0.050 * dx;
	maxx = maxx + 0.050 * dx;
	miny = miny - 0.050 * dy;
	maxy = maxy + 0.050 * dy;

	dx = maxx - minx;
	dy = maxy - miny;

	sprintf( title, "Limits are %+6.4E to %+6.4E on x, %+6.4E to %+6.4E on y", minx, maxx, miny, maxy );	


	if ( FIRST == 0 )
	{
        minsize( MINX_SIZE, MINY_SIZE );
        winid = winopen( title );
        if ( winid == -1 )
                {
                        printf("\033[37m\033[1mCan't open graphics window. Abort.\033[0m\n");
                        myexit( 1 );
		}
		
        winset ( winid );
        FIRST = 1;
        }
        
	ortho2( minx, maxx, miny, maxy);
	doublebuffer();
	gconfig();
	color( 0 );
	clear();
        if ( LARGE_LABELS == 0 )
            loadXfont(4711, "fixed");
        else
            loadXfont(4711, "10x20");
        font(4711);
        color(7);
        cmov2( minx, miny );
/*        charstr( "Data loaded. Rendering ..." );
*/        gflush();
        gsync();
	swapbuffers();

	qdevice(QKEY);
	qdevice(DKEY);
	qdevice(FKEY);
	qdevice(LKEY);
	qdevice(PKEY);
        qdevice(PAGEUPKEY);	
        qdevice(PAGEDOWNKEY);	
	qdevice(DOWNARROWKEY);	
	qdevice(UPARROWKEY);	
	qdevice(RIGHTARROWKEY);	
	qdevice(LEFTARROWKEY);	
	qdevice(LEFTMOUSE);
	qdevice(RIGHTMOUSE);
	do_plot_two();

	return;
}




void 	do_plot_two()
{
	int i;
	Screencoord     left, right, bottom, top;


	if ( have_plot_labels == 0 && PLOT_LABELS == 1 )
	{
	float	stepx, stepy;
	float 	roundminy, roundminx;
	
        gflush();
        gsync();
        backbuffer( 1 );

	have_plot_labels = 1;
	
	color( 6 );
	linewidth( 1 );
	mapcolor( 201  ,50, 50, 50 );
	deflinestyle( 2 , 0x0101 );
	setlinestyle( 2 );

	stepx = pow( 10.0 , floor( log10( (maxx-minx)/5 )) );
	stepy = pow( 10.0 , floor( log10( (maxy-miny)/5 )) );
	
	stepx *= (int)( ( (maxx-minx) / stepx ) / 5 );
	stepy *= (int)( ( (maxy-miny) / stepy ) / 5 );
	
	roundminx = stepx * (int)(minx / stepx);
	roundminy = stepy * (int)(miny / stepy);
	
	for ( i=0 ; i < (int)((maxy-miny)/stepy+2.5)  ; i++ )
	  {
	    if ( i==0 && roundminy < miny ) 
	    		;
	    else
	    {
	    if ( PLOT_LIGHTGRID == 0 )
	    {	    
	    move2( minx, i*stepy+roundminy );
	    draw2( maxx, i*stepy+roundminy );
	    }
	    else
	    {
	    color( 201 );
	    setlinestyle( 0 );
	    move2( minx, i*stepy+roundminy );
	    draw2( maxx, i*stepy+roundminy );
	    setlinestyle( 2 );
	    color( 6 );
	    }
	    
	    if ( i*stepy+roundminy < 0.9850*(maxy-miny)+miny )
	    {
	    cmov2( minx, i*stepy+roundminy );
	    if ( fabs( maxy-miny ) > 9999999 || fabs(maxy) > 9999999 || fabs(miny) > 9999999 )
	    	sprintf( label, "%+6.4E", i*stepy+roundminy );
	    else if ( fabs( maxy-miny ) > 50 )
	    	sprintf( label, "%.0f", i*stepy+roundminy );
	    else if ( fabs( maxy-miny ) > 5 )
	    	sprintf( label, "%.1f", i*stepy+roundminy );
	    else if ( fabs( maxy-miny ) > 1 ) 
	        sprintf( label, "%.2f", i*stepy+roundminy );
	    else
	      	sprintf( label, "%+6.4E", i*stepy+roundminy );
	    charstr( label );
	    
	    if ( AUTOSCALE == YES )
	    {
	    if ( fabs( maxz-minz ) > 9999999 || fabs(maxz) > 9999999 || fabs(minz) > 9999999 )
	    	sprintf( label2, "%+6.4E", ((i*stepy+roundminy)-transsecond)/scalesecond );
	    else if ( fabs( maxy-miny ) > 50 )
	    	sprintf( label2, "%.0f", ((i*stepy+roundminy)-transsecond)/scalesecond );
	    else if ( fabs( maxy-miny ) > 5 )
	    	sprintf( label2, "%.1f", ((i*stepy+roundminy)-transsecond)/scalesecond );
	    else if ( fabs( maxy-miny ) > 1 ) 
	        sprintf( label2, "%.2f", ((i*stepy+roundminy)-transsecond)/scalesecond );
	    else
	      	sprintf( label2, "%+6.4E", ((i*stepy+roundminy)-transsecond)/scalesecond );
            
            getviewport( &left, &right, &bottom, &top );
	    cmov2( maxx-(strwidth(label2)*(maxx-minx)/(right-left)), i*stepy+roundminy );
	    charstr( label2 );
	    }
	    }
	    }
	  }
	for ( i=0 ; i < (int)((maxx-minx)/stepx+2.5) ; i++ )
	  {
	    if ( i==0 && roundminx < minx )
	    		;
	    else
	    {
	    if ( PLOT_LIGHTGRID == 0 )
	    {
	    move2( i*stepx+roundminx, miny );
	    draw2( i*stepx+roundminx, maxy );
	    }
	    else
	    {
	    color( 201 );
	    setlinestyle( 0 );
	    move2( i*stepx+roundminx, miny );
	    draw2( i*stepx+roundminx, maxy );
	    setlinestyle( 2 );
	    color( 6 );
	    }
	    
	    if ( i*stepx+roundminx < 0.9750*(maxx-minx)+minx )
	    {
	    cmov2( i*stepx+roundminx, miny );
	    if ( fabs( maxx-minx ) > 9999999 || fabs(maxx) > 9999999 || fabs(minx) > 9999999 )
	    	sprintf( label, "%+6.4E", i*stepx+roundminx );
	    else if ( fabs( maxx-minx ) > 50 )
	    	sprintf( label, "%.0f", i*stepx+roundminx );
	    else if ( fabs( maxx-minx ) > 5 )
	    	sprintf( label, "%.1f", i*stepx+roundminx );
	    else if ( fabs( maxx-minx ) > 1 ) 
	        sprintf( label, "%.2f", i*stepx+roundminx );
	    else
	      	sprintf( label, "%+6.4E", i*stepx+roundminx );
	    charstr( label );
	    }
	    }
	  }
        setlinestyle( 0 );

	
	color( 7 );
	linewidth( 2 );
	if ( minx < 0 && maxx > 0 )
		{
			move2( 0.0, miny );
			draw2( 0.0, maxy );
		}

	if ( miny < 0 && maxy > 0 )
		{
			move2( minx, 0.0 );
			draw2( maxx, 0.0 );
		}

        gflush();
        gsync();
        swapbuffers();
	return;                			

	}

        frontbuffer( 1 );
        if ( COLOR == YES )
          color(0);
        else
          color(7);
        cmov2( minx, miny );
        charstr( "Rendering ..." );
        gflush();
        gsync();
        backbuffer( 1 );
	
        color( 0 );
        clear();



	if ( PLOT_LABELS == 1 )
	{
	float	stepx, stepy;
	float 	roundminy, roundminx;
	
	have_plot_labels = 1;
	
	color( 6 );
	linewidth( 1 );
	mapcolor( 201  ,50, 50, 50 );
	deflinestyle( 2 , 0x0101 );
	setlinestyle( 2 );

	stepx = pow( 10.0 , floor( log10( (maxx-minx)/5 )) );
	stepy = pow( 10.0 , floor( log10( (maxy-miny)/5 )) );
	
	stepx *= (int)( ( (maxx-minx) / stepx ) / 5 );
	stepy *= (int)( ( (maxy-miny) / stepy ) / 5 );
	
	roundminx = stepx * (int)(minx / stepx);
	roundminy = stepy * (int)(miny / stepy);
	
	for ( i=0 ; i < (int)((maxy-miny)/stepy+2.5)  ; i++ )
	  {
	    if ( i==0 && roundminy < miny ) 
	    		;
	    else
	    {
	    if ( PLOT_LIGHTGRID == 0 )
	    {	    
	    move2( minx, i*stepy+roundminy );
	    draw2( maxx, i*stepy+roundminy );
	    }
	    else
	    {
	    color( 201 );
	    setlinestyle( 0 );
	    move2( minx, i*stepy+roundminy );
	    draw2( maxx, i*stepy+roundminy );
	    setlinestyle( 2 );
	    color( 6 );
	    }
	    
	    if ( i*stepy+roundminy < 0.9850*(maxy-miny)+miny )
	    {
	    cmov2( minx, i*stepy+roundminy );
	    if ( fabs( maxy-miny ) > 9999999 || fabs(maxy) > 9999999 || fabs(miny) > 9999999 )
	    	sprintf( label, "%+6.4E", i*stepy+roundminy );
	    else if ( fabs( maxy-miny ) > 50 )
	    	sprintf( label, "%.0f", i*stepy+roundminy );
	    else if ( fabs( maxy-miny ) > 5 )
	    	sprintf( label, "%.1f", i*stepy+roundminy );
	    else if ( fabs( maxy-miny ) > 1 ) 
	        sprintf( label, "%.2f", i*stepy+roundminy );
	    else
	      	sprintf( label, "%+6.4E", i*stepy+roundminy );
	    charstr( label );
	    }

	    if ( AUTOSCALE == YES )
	    {
	    if ( fabs( maxz-minz ) > 9999999 || fabs(maxz) > 9999999 || fabs(minz) > 9999999 )
	    	sprintf( label2, "%+6.4E", ((i*stepy+roundminy)-transsecond)/scalesecond );
	    else if ( fabs( maxy-miny ) > 50 )
	    	sprintf( label2, "%.0f", ((i*stepy+roundminy)-transsecond)/scalesecond );
	    else if ( fabs( maxy-miny ) > 5 )
	    	sprintf( label2, "%.1f", ((i*stepy+roundminy)-transsecond)/scalesecond );
	    else if ( fabs( maxy-miny ) > 1 ) 
	        sprintf( label2, "%.2f", ((i*stepy+roundminy)-transsecond)/scalesecond );
	    else
	      	sprintf( label2, "%+6.4E", ((i*stepy+roundminy)-transsecond)/scalesecond );
            
            getviewport( &left, &right, &bottom, &top );
	    cmov2( maxx-(strwidth(label2)*(maxx-minx)/(right-left)), i*stepy+roundminy );
	    charstr( label2 );
	    }

	    }
	  }
	for ( i=0 ; i < (int)((maxx-minx)/stepx+2.5) ; i++ )
	  {
	    if ( i==0 && roundminx < minx )
	    		;
	    else
	    {
	    if ( PLOT_LIGHTGRID == 0 )
	    {
	    move2( i*stepx+roundminx, miny );
	    draw2( i*stepx+roundminx, maxy );
	    }
	    else
	    {
	    color( 201 );
	    setlinestyle( 0 );
	    move2( i*stepx+roundminx, miny );
	    draw2( i*stepx+roundminx, maxy );
	    setlinestyle( 2 );
	    color( 6 );
	    }
	    
	    if ( i*stepx+roundminx < 0.9750*(maxx-minx)+minx )
	    {
	    cmov2( i*stepx+roundminx, miny );
	    if ( fabs( maxx-minx ) > 9999999 || fabs(maxx) > 9999999 || fabs(minx) > 9999999 )
	    	sprintf( label, "%+6.4E", i*stepx+roundminx );
	    else if ( fabs( maxx-minx ) > 50 )
	    	sprintf( label, "%.0f", i*stepx+roundminx );
	    else if ( fabs( maxx-minx ) > 5 )
	    	sprintf( label, "%.1f", i*stepx+roundminx );
	    else if ( fabs( maxx-minx ) > 1 ) 
	        sprintf( label, "%.2f", i*stepx+roundminx );
	    else
	      	sprintf( label, "%+6.4E", i*stepx+roundminx );
	    charstr( label );
	    }
	    }
	  }
        setlinestyle( 0 );
        }
        	else
        		have_plot_labels = 0;



	linewidth( 1 );
	color( 3 );


	if ( ERROR_BARS == YES )
	{

	mapcolor( 201  ,150, 150, 150 );	
	color( 201 );
	
	getviewport( &left, &right, &bottom, &top );
	for ( i=0 ; i < N ; i++ )
		{
		move2( x[i], y[i]-z[i] );
		draw2( x[i], y[i]+z[i] );
		move2( x[i]-4.0*dx/(right-left) , y[i]-z[i] );
		draw2( x[i]+4.0*dx/(right-left) , y[i]-z[i] );
		move2( x[i]-4.0*dx/(right-left) , y[i]+z[i] );
		draw2( x[i]+4.0*dx/(right-left) , y[i]+z[i] );
		}
	

        mapcolor( 200  ,243, 146, 46 );
	color( 200 );

	if ( dots == NO && filled_dots == NO )
	{
	move2( x[0], y[0] );
	for ( i=1 ; i < N ; i++ )
		draw2( x[i], y[i] );
	}
	else if ( filled_dots == YES )
	{
	getviewport( &left, &right, &bottom, &top );

	for ( i=0 ; i < N ; i++ )
		arcxf( x[i], y[i], 4.0*dx/(right-left), 4.0*dy/(top-bottom), 0, 3600 );
        
        if ( dots == NO )
        {
	move2( x[0], y[0] );
	for ( i=1 ; i < N ; i++ )
		draw2( x[i], y[i] );
        }                
        
	}
	else
	{
	for ( i=0 ; i < N ; i++ )
		pnt2( x[i], y[i] );
	}

	
	

	}
	else
	{


	if ( dots == NO && filled_dots == NO )
	{

            i = 0;
            while ( x[i] < minx && i < N )
                i++;
       
            if ( i > 0 )
                i--;

	    move2( x[i], z[i] );
	    while ( x[i] < maxx && i < N )
                {
		    draw2( x[i], z[i] );
                    i++;
                }

            if ( i < N )
                draw2( x[i], z[i] );
	}
	else if ( filled_dots == YES )
	{
	getviewport( &left, &right, &bottom, &top );

            i = 0;
            while ( x[i] < minx && i < N )
                i++;
            
            if ( i > 0 )
                i--;

	    while ( x[i] < maxx && i < N )
                {
		    arcxf( x[i], z[i], 4.0*dx/(right-left), 4.0*dy/(top-bottom), 0, 3600 );
                    i++;
                }

            if ( i < N )
                arcxf( x[i], z[i], 4.0*dx/(right-left), 4.0*dy/(top-bottom), 0, 3600 );

        
        
        if ( dots == NO )
        {

                i = 0;
                while ( x[i] < minx && i < N )
                    i++;
        
                if ( i > 0 )
                    i--;

                move2( x[i], z[i] );
                while ( x[i] < maxx && i < N )
                {
                    draw2( x[i], z[i] );
                    i++;
                }

                if ( i < N )
                    draw2( x[i], z[i] );
        }                
        
	}
	else
	{

            i = 0;
            while ( x[i] < minx && i < N )
                i++;
        
            if ( i > 0 )
                i--;

	    while ( x[i] < maxx && i < N )
                {
		    pnt2( x[i], z[i] );
                    i++;
                }

            if ( i < N )
                pnt2( x[i], z[i] );

	}
    


	
        mapcolor( 200  ,243, 146, 46 );
	color( 200 );


	if ( dots == NO && filled_dots == NO )
	{
            i = 0;
            while ( x[i] < minx && i < N )
                i++;
       
            if ( i > 0 )
                i--;

	    move2( x[i], y[i] );
	    while ( x[i] < maxx && i < N )
                {
		    draw2( x[i], y[i] );
                    i++;
                }

            if ( i < N )
                draw2( x[i], y[i] );
	}
	else if ( filled_dots == YES )
	{
	    getviewport( &left, &right, &bottom, &top );

            i = 0;
            while ( x[i] < minx && i < N )
                i++;
            
            if ( i > 0 )
                i--;

	    while ( x[i] < maxx && i < N )
                {
		    arcxf( x[i], y[i], 4.0*dx/(right-left), 4.0*dy/(top-bottom), 0, 3600 );
                    i++;
                }

            if ( i < N )
                arcxf( x[i], y[i], 4.0*dx/(right-left), 4.0*dy/(top-bottom), 0, 3600 );
            
            if ( dots == NO ) 
            {
                i = 0;
                while ( x[i] < minx && i < N )
                    i++;
        
                if ( i > 0 )
                    i--;

                move2( x[i], y[i] );
                while ( x[i] < maxx && i < N )
                {
                    draw2( x[i], y[i] );
                    i++;
                }

                if ( i < N )
                    draw2( x[i], y[i] );
            }
	}
	else
	{

            i = 0;
            while ( x[i] < minx && i < N )
                i++;
        
            if ( i > 0 )
                i--;

	    while ( x[i] < maxx && i < N )
                {
		    pnt2( x[i], y[i] );
                    i++;
                }

            if ( i < N )
                pnt2( x[i], y[i] );
	}
        }

	
	color( 7 );
	linewidth( 2 );
	if ( minx < 0 && maxx > 0 )
		{
			move2( 0.0, miny );
			draw2( 0.0, maxy );
		}

	if ( miny < 0 && maxy > 0 )
		{
			move2( minx, 0.0 );
			draw2( maxx, 0.0 );
		}

        gflush();
        gsync();
	swapbuffers();

	return;                			
}





/********
*
*
* Bicubic interpolation
*
* Given   1. A global 2D matrix with the name 'mat' (first index treated as x)
*         2. The coordinates (in zero-based) grid units of the unit square of interest, and,
*         3. Two values for x & y both within the [0,1] unit square limits
*
* return  the Bicubic estimator p(x,y)
*
*
*
********/

float bicubic( gx, gy, x, y )
int	gx, gy;
double	x, y;
{
  double	an1, a0, a1, a2;
  double	bn1, b0, b1, b2;

  an1 = mat[gx-1][gy-1];
  a0  = mat[gx  ][gy-1];
  a1  = mat[gx+1][gy-1];
  a2  = mat[gx+2][gy-1];
  bn1 = ( 2*a0 + x*(-an1+a1) + x*x*(2*an1-5*a0+4*a1-a2) + x*x*x*(-an1+3*a0-3*a1+a2) ) / 2;

  an1 = mat[gx-1][gy];
  a0  = mat[gx  ][gy];
  a1  = mat[gx+1][gy];
  a2  = mat[gx+2][gy];
  b0  = ( 2*a0 + x*(-an1+a1) + x*x*(2*an1-5*a0+4*a1-a2) + x*x*x*(-an1+3*a0-3*a1+a2) ) / 2;

  an1 = mat[gx-1][gy+1];
  a0  = mat[gx  ][gy+1];
  a1  = mat[gx+1][gy+1];
  a2  = mat[gx+2][gy+1];
  b1  = ( 2*a0 + x*(-an1+a1) + x*x*(2*an1-5*a0+4*a1-a2) + x*x*x*(-an1+3*a0-3*a1+a2) ) / 2 ;

  an1 = mat[gx-1][gy+2];
  a0  = mat[gx  ][gy+2];
  a1  = mat[gx+1][gy+2];
  a2  = mat[gx+2][gy+2];
  b2  = ( 2*a0 + x*(-an1+a1) + x*x*(2*an1-5*a0+4*a1-a2) + x*x*x*(-an1+3*a0-3*a1+a2) ) / 2 ;


  return( ( 2*b0 + y*(-bn1+b1) + y*y*(2*bn1-5*b0+4*b1-b2) + y*y*y*(-bn1+3*b0-3*b1+b2) ) / 2 );
}




void contours()
{
	char	line[500000];	
	char 	*p;
	float	junk, val;
	int	i, k;
	int	ii, kk=0;
	double	xfrac, yfrac;
	FILE 	*out;
	float	scale;
	int	wx, wy;
	short 	data;
	long	dev;




	while( fgets( line, 499999, stdin ) != NULL )
	{
	if ( strlen( line ) > 499997 )
		{
			printf("\033[31m\033[1mLine longer than half a million characters ? Abort.\033[0m\n");
			myexit(1);
		}
		
	if ( (columns = sscanf( line, "%f %f %f %f", &junk, &junk, &junk, &junk)) >= 1 )
		break;
	
	printf("\033[37m\033[1mCaution: header line skipped:\033[0m %s", line );	
	}	

	/* Number of columns */
	p = &line[0];
	columns = 0;
	while ( sscanf( p, "%f%n", &junk, &i) == 1 )
		{
			columns++;
			p += i;
		}

	mat   = matrix( 0, columns+2, 0, MAXML );
	inter = matrix( 0, columns+2, 0, MAXML );


	/* Place the first line in matrix */
	p = &line[0];
	k = 1;
	while ( sscanf( p, "%f%n", &val, &i) == 1 )
		{
			mat[k][1] = val;
			k++;
			p += i;
		}




	/* Read the rest */
	lines = 1;				
	while( 1 )
	{
		for ( i=0 ; i < columns ; i++ )
			{
				if ( scanf("%f", &val) != 1 )
					break;
                                if ( isnan(val) || isinf(val) )
                                  {
                                    printf("\033[31m\033[1mData appear to contain one or more 'NaN' and/or 'Inf'. Goodbye.\033[0m\n");
                                    myexit(1);
                                  }
				mat[i+1][lines+1] = val;
			}

		if ( i == 0 )
			break;
		else if ( i != columns )
			{
				printf("\033[31m\033[1mNumber of columns in matrix not constant ? Abort.\033[0m\n");
				myexit(1);
			}
		lines++;

		if ( lines > MAXML - 2 )
			{
                            mat_inter   = mat;
                            inter_inter = inter;

                            mat   = matrix( 0, columns+2, 0, MAXML+STEPM );
                            inter = matrix( 0, columns+2, 0, MAXML+STEPM );      


                            for ( ii=0 ; ii <= columns ; ii++ )
                            for ( kk=0 ; kk <= lines ; kk++ )
                                {
                                  mat[ii][kk] = mat_inter[ii][kk];                                                        
                                }
                                
                            free_matrix(   mat_inter, 0, columns+2, 0, MAXML );
                            free_matrix( inter_inter, 0, columns+2, 0, MAXML );

                            MAXML += STEPM;
			}
	}

	

        if ( LOG == YES || LOGLOG == YES )
        {	
	for ( i=1 ; i <= columns ; i++ )
	for ( k=1 ; k <= lines   ; k++ )
	  {
	    if ( mat[i][k] <= 0 )
	      {
		printf("\033[31m\033[1mLogarithm requested but negative values present ? Abort.\033[0m\n");
		myexit(1);
	      }
            mat[i][k] = log( mat[i][k] );
          }
        }



	/* Sub-matrix requested ? */
	if ( SUBMATRIX == YES )
	{
	        FILE	*sub;
		int	sourcec, sourcel;
		
		
		for ( ii=1 ; ii <= newc ; ii++ )
		for ( kk=1 ; kk <= newl ; kk++ )
			{
				sourcec = oric + ii -1;
				sourcel = oril + kk -1;
				
				if ( sourcec < 1 || sourcec > columns || sourcel < 1 || sourcel > lines )
					{
					printf("\033[31m\033[1mWrong size-origin combination for submatrix.\033[0m\n");
					myexit(1);
					}
				
				mat[ii][kk] = mat[sourcec][sourcel];
				
			}


		sub = fopen("submatrix.dat", "w" );
                if ( sub != NULL )
                {		
		    for ( ii=1 ; ii <= newc ; ii++ )
		      {
		        for ( kk=1 ; kk <= newl ; kk++ )
		            fprintf(sub, " %+10.8E", mat[ii][kk] );
                
                        fprintf(sub, "\n");
                      }
                    fclose( sub );
		}
		
		columns = newc;
		lines = newl;
	}


	/* Strange matrix ? */
	if ( columns <= 1 || lines <= 1 )
		{
			printf("\033[31m\033[1mNothing to plot.\033[0m\n");
			myexit(1);
		}
        


	if ( ((float)(columns) / (lines) > 10.0 || (float)(columns) / (lines) < 0.10) && SECTIONS < 2 )
		{
		printf("\033[37m\033[1mNice matrix (%dx%d). Expect problems ...\033[0m\n", columns, lines);
		}



	/* Copy first-last lines-columns for bicubic interpolation */
	for ( i=1 ; i <= columns ; i++ )
		{
			mat[i][0] = mat[i][1];
			mat[i][lines+1] = mat[i][lines];
		}
	for ( i=1 ; i <= lines ; i++ )
		{
			mat[0][i] = mat[1][i];
			mat[columns+1][i] = mat[columns][i];
		}
	
	mat[0][0] = mat[1][1];
	mat[0][lines+1] = mat[1][lines];
	mat[columns+1][0] = mat[columns][1];
	mat[columns+1][lines+1] = mat[columns][lines];
	

/* 	for ( i=0 ; i <= lines+1 ; i++ )
		{
			for ( k=0 ; k <= columns+1 ; k++ )
				printf(" %8.4f", mat[k][i]);
			printf("\n");
		}
*/


	/* For tiny matrices, use bicubic interpolation */
        mult = 1;
	ori_col = columns;
	ori_lin = lines;
	if ( columns < SMALLM || (lines/SECTIONS) < SMALLM )
	{
		if ( (SMALLM / (lines/SECTIONS)) > (SMALLM / columns) )
			mult = (int)( SMALLM / (columns-1) + 0.50) + 1;
		else
			mult = (int)( SMALLM / ((lines/SECTIONS)-1) + 0.50) + 1;
		
		if ( VERBOSE )			
		printf("\033[37m\033[1mBicubic interpolation to %dx%d ...\033[0m\n", mult*(columns-1), mult*(lines-1) );


		/* Re-allocate memory for big matrix */
                mat_inter   = mat;
                inter_inter = inter;

                if ( SECTIONS <= 1 )
                {
                mat   = matrix( 0, 2*mult*(columns-1)+1, 0, 2*mult*(lines-1)+1 );
                inter = matrix( 0, 2*mult*(columns-1)+1, 0, 2*mult*(lines-1)+1 );
                }
                else
                {
                mat   = matrix( 0, mult*(columns+1)+1, 0, mult*(lines+1)+1 );
                inter = matrix( 0, mult*(columns+1)+1, 0, mult*(lines+1)+1 );      
                }

                for ( ii=0 ; ii <= columns+1 ; ii++ )
                for ( kk=0 ; kk <= lines+1 ; kk++ )
                    {
                      mat[ii][kk] = mat_inter[ii][kk];                                                        
                    }
                    
                free_matrix(   mat_inter, 0, columns+2, 0, MAXML );
                free_matrix( inter_inter, 0, columns+2, 0, MAXML );


		/* indeces for interpolated matrix */
		for ( ii=0 ; ii < mult*(columns-1) ; ii++ )
		for ( kk=0 ; kk < mult*(lines-1)   ; kk++ )
		{
			xfrac = (double)(ii % mult) / (double)(mult);
			yfrac = (double)(kk % mult) / (double)(mult);
			inter[ii][kk] = bicubic( ii/mult+1, kk/mult+1, xfrac, yfrac );
		}
	
	
		if ( VERBOSE )
		{
		/* Write out */	
		out = fopen("bicubic.dat", "w" );
		if ( out != NULL )
		{
		for ( i=0 ; i < mult*(lines-1) ; i++ )
			{
				for ( k=0 ; k < mult*(columns-1) ; k++ )
					fprintf(out, " %10.4f", inter[k][i]);
				fprintf(out, "\n");
			}
		printf("\033[37m\033[1mInterpolated matrix written in file 'bicubic.dat' ...\033[0m\n");
		}
		fclose ( out );

/*
		out = fopen("bicubic_padded.dat", "w" );
		if ( out != NULL )
		{
		int	mindex;
		
		if ( lines > columns)
			mindex = mult*(lines-1);
		else
			mindex = mult*(columns-1);
		
		for ( i=0 ; i < mindex ; i++ )
			{
				for ( k=0 ; k < mindex ; k++ )
					{
						if ( k < mult*(columns-1) && i < mult*(lines-1) )
							fprintf(out, " %10.4f", inter[k][i]);
						else
							fprintf(out, " %10.4f", inter[0][0]);
					}
				fprintf(out, "\n");
			}
		printf("\033[37m\033[1mPadded (square) matrix written in file 'bicubic_padded.dat' ...\033[0m\n");
		}
		fclose ( out );
*/


		}
		

		/* Update number of columns and lines */
                columns = mult*(columns-1);
		lines = mult*(lines-1);

	}
	else
	{

	        if ( lines > 2*YMAXSCREEN || columns > 2*XMAXSCREEN )
	          BICUB = NO;
                else
                {
	          BICUB = YES;
	        
	          /* Increase matrices sizes for zooming ... */
	          mat_inter   = mat;
	          inter_inter = inter;

	          mat   = matrix( 0, 2*(columns+2), 0, 2*(lines+2) );
	          inter = matrix( 0, 2*(columns+2), 0, 2*(lines+2) );      


	          for ( ii=0 ; ii <= columns+1 ; ii++ )
	          for ( kk=0 ; kk <= lines+1 ; kk++ )
                      {
                        mat[ii][kk] = mat_inter[ii][kk];                                                        
                      }
                    
                  free_matrix(   mat_inter, 0, columns+2, 0, MAXML );
                  free_matrix( inter_inter, 0, columns+2, 0, MAXML );


                  mult *= 2;
                  for ( ii=0 ; ii < mult*(ori_col-1) ; ii++ )
                  for ( kk=0 ; kk < mult*(ori_lin-1)   ; kk++ )
                      {
		  	xfrac = (double)(ii % mult) / (double)(mult);
		  	yfrac = (double)(kk % mult) / (double)(mult);
			inter[ii][kk] = bicubic( ii/mult+1, kk/mult+1, xfrac, yfrac );
                          
                      }

                  mult /= 2;
                  for ( ii=0 ; ii < mult*(ori_col-1) ; ii++ )
                  for ( kk=0 ; kk < mult*(ori_lin-1)   ; kk++ )
                    {
			xfrac = (double)(ii % mult) / (double)(mult);
			yfrac = (double)(kk % mult) / (double)(mult);
			mat[ii][kk] = bicubic( ii/mult+1, kk/mult+1, xfrac, yfrac );
                    }
		}


                for ( ii=0 ; ii < columns ; ii++ )
                for ( kk=0 ; kk < lines   ; kk++ )
		  inter[ii][kk] = mat[ii+1][kk+1];
	}


	
	/* OK. Have the possibly interpolated matrix in inter[][] and dimensions. Get on with it ... */

	/* Get min, max, mean and rmsd */	
	min = inter[0][0];
	max = inter[0][0];
	mean = 0.0;
	rmsd = 0.0;
	for ( ii=0 ; ii < columns ; ii++ )
	for ( kk=0 ; kk < lines   ; kk++ )
		{
		        val = inter[ii][kk];
			if ( val > max )
				max = val;
			if ( val < min )
				min = val;
			mean += val;
		}	
	mean /= (ii*kk);

	

	for ( ii=0 ; ii < columns ; ii++ )
	for ( kk=0 ; kk < lines   ; kk++ )
		{
                        val = inter[ii][kk];
			rmsd += (val -mean)*(val-mean);
		}


	rmsd = sqrt( rmsd / (ii*kk));

	if ( HAVE_MIN_MAX == YES && DRAW_HISTOGRAM == NO )
	  {
	    min = MIN_FROM_CML;
	    max = MAX_FROM_CML;
	  }


	c_step = 1.0;

	if ( VERBOSE)
		{
		printf("\033[37m\033[1mMin, max, mean, rmsd are %f %f %f %f\033[0m\n", min, max, mean, rmsd);
		printf("\033[37m\033[1mFirst contour at %f and then every %f\033[0m\n", mean+rmsd/2, c_step * rmsd/2);
		}


        lines = (lines - (ori_lin%SECTIONS>0?0:1)*(SECTIONS-1)*mult) / SECTIONS;


	/* Initial window size in pixels ... */
	scale = 1.0;
	if ( columns > lines )
	{
		if ( columns > MINS_SIZE )
			scale = (float)(MINS_SIZE) / columns;
	}
	else
	{
		if ( lines > MINS_SIZE )
			scale = (float)(MINS_SIZE) / lines;
	}	

	wx = (int)(scale * columns + 0.50);
	wy = (int)(scale * lines   + 0.50);
	
	sprintf( title, "Min: %+6.4E Max: %+6.4E Mean: %+6.4E Sigma: %+6.4E", min, max, mean, rmsd );

        minsize( wx, wy );

	if ( COLOR == YES )
		maxsize( (int)(columns), (int)(lines) );

	keepaspect( wx, wy );
		
        winid = winopen( title );
        if ( winid == -1 )
                {
                        printf("\033[31m\033[1mCan't open graphics window. Abort.\033[0m\n");
                        myexit( 1 );
		}
		
        winset ( winid );
        deflinestyle( 1, 257);

        if ( COLOR == YES )
        {
        	int	fraction;
        	int	r, g, b;

        	for ( fraction=1 ; fraction <= 255 ; fraction++ )
        		{
        			getcol( (float)(fraction) / 255.0, &r, &g, &b );
        			mapcolor( fraction, r, g, b );
			}        		
        }
        

        dx = columns;
        dy = lines;
	ortho2( 0, columns, lines, 0);
	minx = 0;
	miny = lines;
	doublebuffer();
	gconfig();
        color( 0 );
        clear();
	qdevice(QKEY);
	qdevice(PKEY);
	qdevice(LEFTMOUSE);
	qdevice(NKEY);
        qdevice(ZKEY);
	qdevice(CKEY);
	qdevice(MINUSKEY);
	qdevice(EQUALKEY);
        qdevice(DOWNARROWKEY);
        qdevice(UPARROWKEY);
        qdevice(PAGEUPKEY);
        qdevice(PAGEDOWNKEY);


        if ( LARGE_LABELS == 0 )
            loadXfont(4711, "fixed");
        else
            loadXfont(4711, "10x20");
        font(4711);
        if ( COLOR == YES )
          color(120);
        else
          color(7);
        cmov2( 0, lines );
/*        charstr( "Matrix loaded. Rendering ..." );
*/        gflush();
        gsync();
	swapbuffers();


	do_plot_contours();

	while ( (dev = qread(&data) ) ) 
		{
			if ( dev == REDRAW )
				{
					if ( FIRST_REDRAW == NO )
					{
					reshapeviewport();
					do_plot_contours();
					}
					else
						FIRST_REDRAW = NO;
				}
			if ( dev == EQUALKEY )
				{
					c_step /= 2;
					reshapeviewport();
					do_plot_contours();
					usleep(500000);
					qreset();
				}
			if ( dev == MINUSKEY )
				{
					c_step *= 2;
					reshapeviewport();
					do_plot_contours();
					usleep(500000);
					qreset();
				}
                        if ( dev == DOWNARROWKEY )
                                {
					SECTION++;
                                        if ( SECTION >= SECTIONS )
                                            SECTION = SECTIONS-1;
					reshapeviewport();
					do_plot_contours();
                                        usleep(500000);
					qreset();
                                }
                        if ( dev == UPARROWKEY )
                                {
					SECTION--;
                                        if ( SECTION <= 0 )
                                            SECTION = 0;
					reshapeviewport();
					do_plot_contours();
                                        usleep(500000);
					qreset();
                                }
                        if ( dev == PAGEDOWNKEY )
                                {
					SECTION = SECTIONS-1;;
					reshapeviewport();
					do_plot_contours();
                                        usleep(500000);
					qreset();
                                }
                        if ( dev == PAGEUPKEY )
                                {
					SECTION = 0;
					reshapeviewport();
					do_plot_contours();
                                        usleep(500000);
					qreset();
                                }
			if ( dev == QKEY )
				myexit( 0 );
			if ( dev == PKEY )
				gl2ppm("plot.ppm");
			if ( dev == LEFTMOUSE && IN_ZOOM == NO )
				{
					Screencoord	left, right, bottom, top;
					float	xpos, ypos;
					Icoord	orix, oriy;
					
					xpos = getvaluator(MOUSEX);
					ypos = getvaluator(MOUSEY);
					getorigin( &orix, &oriy );
					getviewport( &left, &right, &bottom, &top );
					
					printf("Fractional : %9.7f %9.7f   Value : %+10.8E\n", 
						(columns*(xpos-orix)/(right-left))/columns, 
						(lines - lines*(ypos-oriy)/(top-bottom))/lines, 
						inter[(int)(columns*(xpos-orix)/(right-left) + 0.50)][(int)((lines - lines*(ypos-oriy)/(top-bottom)) +0.50)]  );
					usleep(500000);
					qreset();
				}
			if ( dev == NKEY )
				{
					if ( DRAW_NEG == NO )
						DRAW_NEG = YES;
					else
						DRAW_NEG = NO;
						
					reshapeviewport();
					do_plot_contours();
					usleep(500000);
					qreset();
				}
			if ( dev == CKEY && COLOR == YES)
				{
					if ( DRAW_CONT == NO )
						DRAW_CONT = YES;
					else
						DRAW_CONT = NO;
						
					reshapeviewport();
					do_plot_contours();
					usleep(500000);
					qreset();
				}
			if ( dev == ZKEY && SECTIONS <= 1 )
				{
					Screencoord	left, right, bottom, top;
					float	xpos, ypos;
					Icoord	orix, oriy;
					float	dx, dy;
	

					if ( IN_ZOOM == NO )
					{
					IN_ZOOM = YES;
        

                                        frontbuffer( 1 );
                                        if ( COLOR == YES )
                                          color(120);
                                        else
                                          color(7);
                                        cmov2( minx, miny );
                                        charstr( "Rendering ..." );
                                        gflush();
                                        gsync();
                                        backbuffer( 1 );
										
					xpos = getvaluator(MOUSEX);
					ypos = getvaluator(MOUSEY);
					getorigin( &orix, &oriy );
					getviewport( &left, &right, &bottom, &top );


                                        if ( BICUB == YES )
                                        {
					dx = 2 * columns*(xpos-orix)/(right-left);
					dy = 2 * (lines - lines*(ypos-oriy)/(top-bottom));

					ortho2( (int)(dx-columns/2), (int)(dx+columns/2), 
						(int)(dy+lines/2), (int)(dy-lines/2));
		   
		                        minx = (int)(dx-columns/2);
		                        miny = (int)(dy+lines/2);

		    			mult *= 2;
		                        
                                        max = bicubic( 1, 1, 1.0, 1.0 );
                                        min = bicubic( 1, 1, 1.0, 1.0 );
					for ( ii=0 ; ii < mult*(ori_col-1) ; ii++ )
					for ( kk=0 ; kk < mult*(ori_lin-1)   ; kk++ )
					{
						xfrac = (double)(ii % mult) / (double)(mult);
						yfrac = (double)(kk % mult) / (double)(mult);
                                                val = bicubic( ii/mult+1, kk/mult+1, xfrac, yfrac );
                                                inter[ii][kk] = val;
                                                if ( val > max )
                                                  max = val;
                                                if ( val < min )
                                                  min = val;
					}
					
					if ( HAVE_MIN_MAX == YES && DRAW_HISTOGRAM == NO )
					{
					  min = MIN_FROM_CML;
                                          max = MAX_FROM_CML;
                                        }
					
					columns = mult*(ori_col-1);
					lines = mult*(ori_lin-1);


                                        }
                                        else
                                        {
					dx = columns*(xpos-orix)/(right-left);
					dy = (lines - lines*(ypos-oriy)/(top-bottom));

					ortho2( (int)(dx-columns/4), (int)(dx+columns/4), 
						(int)(dy+lines/4), (int)(dy-lines/4));
		   
		                        minx = (int)(dx-columns/4);
		                        miny = (int)(dy+lines/4);
                                        }
                                        

                                        
                                        

					reshapeviewport();						
					do_plot_contours();
					usleep(500000);
					qreset();
					}
					else
					{
					IN_ZOOM = NO;

                                        frontbuffer( 1 );
                                        if ( COLOR == YES )
                                          color(120);
                                        else
                                          color(7);
                                        cmov2( minx, miny );
                                        charstr( "Rendering ..." );
                                        gflush();
                                        gsync();
                                        backbuffer( 1 );
					
					if ( BICUB == YES )
					{					

					mult /= 2;
                                        max = bicubic( 1, 1, 1.0, 1.0 );
                                        min = bicubic( 1, 1, 1.0, 1.0 );
					for ( ii=0 ; ii < mult*(ori_col-1) ; ii++ )
					for ( kk=0 ; kk < mult*(ori_lin-1)   ; kk++ )
					{
						xfrac = (double)(ii % mult) / (double)(mult);
						yfrac = (double)(kk % mult) / (double)(mult);
                                                val = bicubic( ii/mult+1, kk/mult+1, xfrac, yfrac );
                                                inter[ii][kk] = val;
                                                if ( val > max )
                                                  max = val;
                                                if ( val < min )
                                                  min = val;
					}

					if ( HAVE_MIN_MAX == YES && DRAW_HISTOGRAM == NO )
					{
					  min = MIN_FROM_CML;
                                          max = MAX_FROM_CML;
                                        }
					
					columns = mult*(ori_col-1);
					lines = mult*(ori_lin-1);

                                        }
                                        

					ortho2( 0, columns, lines, 0 );

                                        minx = 0;
                                        miny = lines;

					reshapeviewport();						
					do_plot_contours();
					usleep(500000);
					qreset();
					
					}
					
				}
			
		}





	myexit(1);
}



void 	do_plot_contours()
{
	int	ii, kk;
	int	v1, v2;
	float	vv1, vv2, prod;

        frontbuffer( 1 );
        if ( COLOR == YES )
          color(120);
        else
          color(7);
        if ( LARGE_LABELS == 0 )
             loadXfont(4711, "fixed");
        else
             loadXfont(4711, "10x20");
        font(4711);
        cmov2( minx, miny );
        if ( SECTIONS > 1 )
            {
                sprintf( label, "Rendering Section %d", SECTION+1 );
                charstr( label );
            }
        else
            charstr( "Rendering ..." );
        gflush();
        gsync();
        backbuffer( 1 );

        color( 0 );
        clear();


        if ( COLOR == YES )
        {
	for ( ii=0 ; ii < columns ; ii++ )
	for ( kk=0 ; kk < lines ; kk++ )
		{	
		        v1 = (int)   ( 255.0 * ( inter[ii][kk+(lines+(ori_lin%SECTIONS>0?0:1)*mult)*SECTION] - min ) / ( max - min ) + 0.50 );
		        if ( v1 > 255 )
		          v1 = 255;
                        if ( v1 < 0 )
                          v1 = 0;
			color( v1 );
			pnt2( ii, kk );
			
		}

	if ( DRAW_CONT == YES )
	{
	color( 0 );
	setlinestyle( 0 );
	linewidth( 1 );
	move2( 0, lines / 2);
	draw2( columns, lines / 2);
	move2( columns/2, 0);
	draw2( columns/2, lines);

	setlinestyle( 1 );
	move2( 0, lines / 4);
	draw2( columns, lines / 4);
	move2( columns/4, 0);
	draw2( columns/4, lines);
	move2( 0, 3 * lines / 4);
	draw2( columns, 3 * lines / 4);
	move2( 3 * columns/4, 0);
	draw2( 3 * columns/4, lines);
	


	for ( ii=0 ; ii < columns-1 ; ii++ )
	for ( kk=0 ; kk < lines   ; kk++ )
		{
			vv1 = inter[ii][kk+(lines+(ori_lin%SECTIONS>0?0:1)*mult)*SECTION];
			vv2 = inter[ii+1][kk+(lines+(ori_lin%SECTIONS>0?0:1)*mult)*SECTION];
			prod = (vv1 - mean) * (vv2 - mean);
			v1 = (int)((vv1 - mean) / (c_step*rmsd/2));
			v2 = (int)((vv2 - mean) / (c_step*rmsd/2));
			
			if ( v1 != v2 && v1 >= 0 && v2 >= 0 )
				pnt2( ii, kk );
			if ( DRAW_NEG == YES && prod < 0 )
				{
					pnt2( ii-1, kk-1 );
					pnt2( ii-1, kk+0 );
					pnt2( ii-1, kk+1 );
					pnt2( ii+0, kk-1 );
					pnt2( ii+0, kk+1 );
					pnt2( ii+1, kk-1 );
					pnt2( ii+1, kk+0 );
					pnt2( ii+1, kk+1 );
				}
			if ( DRAW_NEG == YES && v1 != v2 && v1 <= 0 && v2 <= 0 && getrand < 0.3333333 )
				pnt2( ii, kk );
		}

	for ( ii=0 ; ii < columns   ; ii++ )
	for ( kk=0 ; kk < lines-1   ; kk++ )
		{
			vv1 = inter[ii][kk+(lines+(ori_lin%SECTIONS>0?0:1)*mult)*SECTION];
			vv2 = inter[ii][kk+(lines+(ori_lin%SECTIONS>0?0:1)*mult)*SECTION+1];
			prod = (vv1 - mean) * (vv2 - mean);
			v1 = (int)((vv1 - mean) / (c_step*rmsd/2));
			v2 = (int)((vv2 - mean) / (c_step*rmsd/2));
			
			if ( v1 != v2 && v1 >= 0 && v2 >= 0)
				pnt2( ii, kk );
			if ( DRAW_NEG == YES && prod < 0 )
				{
					pnt2( ii-1, kk-1 );
					pnt2( ii-1, kk+0 );
					pnt2( ii-1, kk+1 );
					pnt2( ii+0, kk-1 );
					pnt2( ii+0, kk+1 );
					pnt2( ii+1, kk-1 );
					pnt2( ii+1, kk+0 );
					pnt2( ii+1, kk+1 );
				}
			if ( DRAW_NEG == YES && v1 != v2 && v1 <= 0 && v2 <= 0 && getrand < 0.3333333 )
				pnt2( ii, kk );

		}

	}

        gflush();
        gsync();
	swapbuffers();
        return;
        }
        
	color( 7 );
	setlinestyle( 0 );
	linewidth( 1 );
	move2( 0, lines / 2);
	draw2( columns, lines / 2);
	move2( columns/2, 0);
	draw2( columns/2, lines);

	setlinestyle( 1 );
	move2( 0, lines / 4);
	draw2( columns, lines / 4);
	move2( columns/4, 0);
	draw2( columns/4, lines);
	move2( 0, 3 * lines / 4);
	draw2( columns, 3 * lines / 4);
	move2( 3 * columns/4, 0);
	draw2( 3 * columns/4, lines);
	


        mapcolor( 200  ,243, 146, 46 );
	color( 200 );


	for ( ii=0 ; ii < columns-1 ; ii++ )
	for ( kk=0 ; kk < lines   ; kk++ )
		{
			vv1 = inter[ii][kk+(lines+(ori_lin%SECTIONS>0?0:1)*mult)*SECTION];
			vv2 = inter[ii+1][kk+(lines+(ori_lin%SECTIONS>0?0:1)*mult)*SECTION];
			prod = (vv1 - mean) * (vv2 - mean);
			v1 = (int)((vv1 - mean) / (c_step*rmsd/2));
			v2 = (int)((vv2 - mean) / (c_step*rmsd/2));
			
			if ( v1 != v2 && v1 >= 0 && v2 >= 0 )
				pnt2( ii, kk );
			if ( DRAW_NEG == YES && prod < 0 )
				{
					pnt2( ii-1, kk-1 );
					pnt2( ii-1, kk+0 );
					pnt2( ii-1, kk+1 );
					pnt2( ii+0, kk-1 );
					pnt2( ii+0, kk+1 );
					pnt2( ii+1, kk-1 );
					pnt2( ii+1, kk+0 );
					pnt2( ii+1, kk+1 );
				}
			if ( DRAW_NEG == YES && v1 != v2 && v1 <= 0 && v2 <= 0 && getrand < 0.3333333 )
				pnt2( ii, kk );
		}

	for ( ii=0 ; ii < columns   ; ii++ )
	for ( kk=0 ; kk < lines-1   ; kk++ )
		{
			vv1 = inter[ii][kk+(lines+(ori_lin%SECTIONS>0?0:1)*mult)*SECTION];
			vv2 = inter[ii][kk+(lines+(ori_lin%SECTIONS>0?0:1)*mult)*SECTION+1];
			prod = (vv1 - mean) * (vv2 - mean);
			v1 = (int)((vv1 - mean) / (c_step*rmsd/2));
			v2 = (int)((vv2 - mean) / (c_step*rmsd/2));
			
			if ( v1 != v2 && v1 >= 0 && v2 >= 0)
				pnt2( ii, kk );
			if ( DRAW_NEG == YES && prod < 0 )
				{
					pnt2( ii-1, kk-1 );
					pnt2( ii-1, kk+0 );
					pnt2( ii-1, kk+1 );
					pnt2( ii+0, kk-1 );
					pnt2( ii+0, kk+1 );
					pnt2( ii+1, kk-1 );
					pnt2( ii+1, kk+0 );
					pnt2( ii+1, kk+1 );
				}
			if ( DRAW_NEG == YES && v1 != v2 && v1 <= 0 && v2 <= 0 && getrand < 0.3333333 )
				pnt2( ii, kk );

		}
			

        gflush();
        gsync();
	swapbuffers();

	return;
}



void getcol( float val, int *r, int *g, int *b )
{

	if ( val <= 0.25 )
		{
			*r = 255 * ( 4.0 * val * (RED_2 - RED_1) + RED_1 ); 
			*g = 255 * ( 4.0 * val * (GREEN_2 - GREEN_1) + GREEN_1 ); 
			*b = 255 * ( 4.0 * val * (BLUE_2 - BLUE_1) + BLUE_1 ); 
		}
	else if ( val <= 0.50 )
		{
			*r = 255 * ( 4.0 * (val - 0.250) * (RED_3 - RED_2) + RED_2 ); 
			*g = 255 * ( 4.0 * (val - 0.250) * (GREEN_3 - GREEN_2) + GREEN_2 ); 
			*b = 255 * ( 4.0 * (val - 0.250) * (BLUE_3 - BLUE_2) + BLUE_2 ); 
		}
	else if ( val <= 0.75 )
		{
			*r = 255 * ( 4.0 * (val - 0.500) * (RED_4 - RED_3) + RED_3 ); 
			*g = 255 * ( 4.0 * (val - 0.500) * (GREEN_4 - GREEN_3) + GREEN_3 ); 
			*b = 255 * ( 4.0 * (val - 0.500) * (BLUE_4 - BLUE_3) + BLUE_3 ); 
		}
	else
		{
			*r = 255 * ( 4.0 * (val - 0.750) * (RED_5 - RED_4) + RED_4 ); 
			*g = 255 * ( 4.0 * (val - 0.750) * (GREEN_5 - GREEN_4) + GREEN_4 ); 
			*b = 255 * ( 4.0 * (val - 0.750) * (BLUE_5 - BLUE_4) + BLUE_4 ); 
		}

        if ( *r < 0 )
          *r = 0;
        if ( *r > 255 )
          *r = 255;
        if ( *g < 0 )
          *g = 0;
        if ( *g > 255 )
          *g = 255;
        if ( *b < 0 )
          *b = 0;
        if ( *b > 255 )
          *b = 255;

}        






/*****************************************************************************************************************************/
/**   Public domain  memory allocation routines from Numerical recipes.                                                     **/
/*****************************************************************************************************************************/

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
        float *v;

        v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
        if (!v) nrerror("allocation failure in vector()");
        return v-nl+NR_END;
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}

void    	nrerror( char error_text[])
{
	printf("\n\n%s\n\n", &error_text[0] );
	myexit( 1 );
}

float	**matrix( long nrl, long nrh, long ncl, long nch )
{
	long 	i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
	float	**m;

	m = (float **)malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;
	m[nrl] = (float *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) nrerror("allocate failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for (i=nrl+1 ; i <= nrh ; i++ )
		m[i] = m[i-1] + ncol;

	return m;
}

float	***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh )
{
	long 	i, j, nrow=nrh-nrl+1, ncol=nch-ncl+1, ndep = ndh-ndl+1;
	float	***t;

	t=(float ***)malloc((size_t)((nrow+NR_END)*sizeof(float**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	t[nrl] = (float **)malloc((size_t)((nrow*ncol+NR_END)*sizeof(float *)));
	if ( !t[nrl] ) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;
	
	t[nrl][ncl] = (float *)malloc((size_t)((nrow*ncol*ndep + NR_END)*sizeof(float)));
	if (!t[nrl][ncl] ) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for( j=ncl+1 ; j <= nch ; j++ )
		t[nrl][j] = t[nrl][j-1] + ndep;
	for ( i=nrl+1 ; i <= nrh ; i++ )
		{
		t[i] = t[i-1] + ncol;
		t[i][ncl] = t[i-1][ncl] + ncol * ndep;
		for ( j=ncl+1 ; j <= nch ; j++ )
			t[i][j] = t[i][j-1] + ndep;
                                              	}

	return t;
}


unsigned short int	***int3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh )
{
	long 	i, j, nrow=nrh-nrl+1, ncol=nch-ncl+1, ndep = ndh-ndl+1;
	unsigned short int	***t;

	t=(unsigned short int ***)calloc((size_t)((nrow+NR_END)), sizeof(unsigned short int**));
	if (!t) nrerror("allocation failure 1 in intf3tensor()");
	t += NR_END;
	t -= nrl;

	t[nrl] = (unsigned short int **)calloc((size_t)((nrow*ncol+NR_END)), sizeof(unsigned short int *));
	if ( !t[nrl] ) nrerror("allocation failure 2 in intf3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;
	
	t[nrl][ncl] = (unsigned short int *)calloc((size_t)((nrow*ncol*ndep + NR_END)), sizeof(unsigned short int));
	if (!t[nrl][ncl] ) nrerror("allocation failure 3 in intf3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for( j=ncl+1 ; j <= nch ; j++ )
		t[nrl][j] = t[nrl][j-1] + ndep;
	for ( i=nrl+1 ; i <= nrh ; i++ )
		{
		t[i] = t[i-1] + ncol;
		t[i][ncl] = t[i-1][ncl] + ncol * ndep;
		for ( j=ncl+1 ; j <= nch ; j++ )
			t[i][j] = t[i][j-1] + ndep;
                                              	}

	return t;
}





void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
        free((FREE_ARG) (m[nrl]+ncl-NR_END));
        free((FREE_ARG) (m+nrl-NR_END));
}

void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
        long ndl, long ndh)
/* free a float f3tensor allocated by f3tensor() */
{
        free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
        free((FREE_ARG) (t[nrl]+ncl-NR_END));
        free((FREE_ARG) (t+nrl-NR_END));
}

void free_int3tensor(int ***t, long nrl, long nrh, long ncl, long nch,
        long ndl, long ndh)
/* free a float f3tensor allocated by f3tensor() */
{
        free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
        free((FREE_ARG) (t[nrl]+ncl-NR_END));
        free((FREE_ARG) (t+nrl-NR_END));
}



/*****************************************************************************************************************************/
/**   END of Public domain  memory allocation routines from Numerical recipes.                                              **/
/*****************************************************************************************************************************/






/*****************************************************************************************************************************/
/**   KZ filter with three iterations						                                            **/
/*****************************************************************************************************************************/

void kzfilter()
{
  int	ii,kk;
  float	sum;
  int	total;
  float	*iter;
  double max, min;
  float	aim;
  char	titlenew[300];
    

  if ( HAVE_DONE_KZ == YES )
    return;

  frontbuffer( 1 );
  if ( COLOR == YES )
    color(0);
  else
    color(1);
  cmov2( minx, miny );
  charstr( "Kolmogorov-Zurbenko filtering ..." );
  gflush();
  gsync();


  kz   = vector( 0, N);
  iter = vector( 0, N);

  memcpy( iter, y, (N+1)*sizeof(float) );

  welford( N, iter );
/*
  printf("\n\nOriginal : %f %f\n\n", welford_mean, welford_variance );
*/
  aim = welford_variance  ;
  

  max = y[0];
  min = y[0];
		
  for ( i=0 ; i < N ; i++ )
	{
		if ( y[i] > max )
			max = y[i];
		if ( y[i] < min )
			min = y[i];
	}

  sort( N, iter-1 );
  halfwidth = (int)( (max-min) / (2*(iter[ (int)(3.0*N/4.0 + 0.5)]-iter[ (int)(N/4.0 +0.5) ]) / pow( N, 1.0l/3.0l)) + 0.50 );

  halfwidth = (int)((float)(N) / halfwidth + 0.50);

  if ( halfwidth < 2 )
  	halfwidth = N / 100 + 1;

  if ( halfwidth < 2 )
  	{
  		halfwidth = 2;
  	}

/*
  halfwidth = N / 100 + 1;
  if ( halfwidth > 1000 )
    halfwidth = 1000;
*/
  
  if ( halfwidth + HalfWidthOffset > 0 && halfwidth + HalfWidthOffset < N )
    halfwidth += HalfWidthOffset;



  total = 0;
  sum = 0.0;
  /* first iteration */
  for ( ii=0 ; ii < N ; ii++ )
  {
    if ( ii-halfwidth < 1 || ii+halfwidth > N-1 )
    {
    total = 0;
    sum = 0.0;
    for ( kk=(ii-halfwidth>=0?ii-halfwidth:0) ; kk <= (ii+halfwidth<N?ii+halfwidth:N-1) ; kk++)
        {
          sum += y[kk];
          total++;
        }
    kz[ii] = sum / total;
    }
    else
    {
    	sum = sum - y[ii-halfwidth-1] + y[ii+halfwidth];
    	kz[ii] = sum / total;
    }
  }

  /* second iteration */
  total = 0;
  sum = 0.0;
  for ( ii=0 ; ii < N ; ii++ )
  {
    if ( ii-halfwidth < 1 || ii+halfwidth > N-1 )
    {
    total = 0;
    sum = 0.0;
    for ( kk=(ii-halfwidth>=0?ii-halfwidth:0) ; kk <= (ii+halfwidth<N?ii+halfwidth:N-1) ; kk++)
        {
          sum += kz[kk];
          total++;
        }
    iter[ii] = sum / total;
    }
    else
    {
    	sum = sum - kz[ii-halfwidth-1] + kz[ii+halfwidth];
    	iter[ii] = sum / total;
    }
  }


  /* third iteration */
  total = 0;
  sum = 0.0;
  for ( ii=0 ; ii < N ; ii++ )
  {
    if ( ii-halfwidth < 1 || ii+halfwidth > N-1 )
    {
    total = 0;
    sum = 0.0;
    for ( kk=(ii-halfwidth>=0?ii-halfwidth:0) ; kk <= (ii+halfwidth<N?ii+halfwidth:N-1) ; kk++)
        {
          sum += iter[kk];
          total++;
        }
    kz[ii] = sum / total;
    }
    else
    {
    	sum = sum - iter[ii-halfwidth-1] + iter[ii+halfwidth];
    	kz[ii] = sum / total;
    }
  }

  free_vector( iter, 0, N );  

  welford( N, kz );

  sprintf(titlenew, "KZ, width %d, variance explained %.2f%%\n", halfwidth, 100.0*welford_variance/aim );
  wintitle( titlenew );

  if ( welford_variance/aim < 0.70 && KZ_ROUNDS < 10 )
  {
  	if ( KZ_ROUNDS == 0 )
  		{
  		KZ_ROUNDS_DS = welford_variance/aim;
	  	KZ_ROUNDS++;
	  	HalfWidthOffset -= halfwidth/2 ;
	  	kzfilter();	
  		}	
  	else
  		{
  			if ( (welford_variance/aim - KZ_ROUNDS_DS) > 0.10 )
  				{
  				KZ_ROUNDS_DS = welford_variance/aim;
			  	KZ_ROUNDS++;
			  	HalfWidthOffset -= halfwidth/2 ;
			  	kzfilter();	
  				}
			else
				{
				KZ_ROUNDS = 10;
				if ( (welford_variance/aim - KZ_ROUNDS_DS) < 0.030 && halfwidth > 1 )
				{
					HalfWidthOffset += halfwidth/2 ;
					kzfilter();
				}
				}
		}
  }
  
  



  cmov2( minx, miny );
  color(0);
  charstr( "Kolmogorov-Zurbenko filtering ..." );
  gflush();
  gsync();

  HAVE_DONE_KZ = YES;
  return;
  
}






/* SORT.C: from "Numerical Recipes in C" p. 247
 * Sorts an array ra[1..n] into ascending numerical order using the Heapsort
 * algorithm. n is input; ra is replaced on output by its sorted
 * rearrangement.
*/
void sort(int n, float *ra)
{
	int l, j, ir, i;
	float rra;

	l = (n >> 1) + 1;
	ir = n;
/* The index l will be decremented from its initial value down to 1 during
 * the "hiring" (heap creation) phase.  Once it reaches 1, the index ir will
 * be decremented from its initial value down to 1 during the "retirement-
 * and-promotion" (heap selection) phase.
*/
	for ( ; ; )
	{
		if (l > 1)					/* still in hiring phase */
			rra = ra[--l];
		else						/* in retirement-and-promotion phase */
		{
			rra = ra[ir];           /* clear a space at end of array */
			ra[ir]=ra[1];			/* retire the top of the heap into it */
			if (--ir == 1) 			/* done with last promotion */
			{
				ra[1] = rra;
				return;
			}						/* end if */
		}							/* end else */
		i = l;						/* whether we are in the hiring phase */
		j = l << 1;					/* or promotion phase, we here set up */
		while ( j <= ir )
		{
			if ( (j < ir) && (ra[j] < ra[j + 1]) )
				++j;				/* compare to the better underling */
			if ( rra < ra[j] )	/* demote rra */
				{
					ra[i] = ra[j];
					j += (i = j);
				}
			else
			    	j = ir + 1;		/* this is rra's level; set j to */
		}                           /* terminate the sift-down */
		ra[i] = rra;				/* put rra into its slot */
	}
}




void welford(int N, float *data)
{
	double	M, S, oldM, x;
	int	i;
	
	M = 0.0l;
	S = 0.0l;
	for ( i=1 ; i <= N ; i++ )
	{
		x = data[i-1];
		oldM = M;
		M += (x-M)/i;
		S += (x-M)*(x-oldM);
	}

	welford_mean = M;
	welford_variance = S/(N-1.0l);
}









/*
**   Discrete fourier transform of real data. NOT FFT !!!!
**   (but with a bit of openmp to make it somewhat faster.
*/

void DFT(int m,float *y1)
{
   long   i,k;
   double arg;
   double cosarg,sinarg;
   float  *x2=NULL,*y2=NULL;
   float  *x1=NULL;
   int	  done = 0;
   
   if ( HAVE_DFT != 0 )
   {
     x2 = calloc(m, sizeof(float));
     memcpy( x2, data_DFT1, m*sizeof(float) );
     memcpy( data_DFT1, data_DFT2, m*sizeof(float) );
     memcpy( data_DFT2, x2, m*sizeof(float) );
     free(x2);
     if ( HAVE_DFT == 1 )
       HAVE_DFT = 2;
     else
       HAVE_DFT = 1;
     return;
   }   
   
   
   x1 = calloc(m, sizeof(float));
   x2 = calloc(m, sizeof(float));
   y2 = calloc(m, sizeof(float));
   if (x1 == NULL || x2 == NULL || y2 == NULL)
                {
                        fprintf(stderr, "\033[31m\033[1mMemory allocation failure in FT.\033[0m\n");
                        HAVE_DFT = 0;
                        return;
		}
   if ( m > 20000 )
   {
     fprintf(stderr, "\033[32m\033[1mPlot is using naive DFT, and not FFT.\033[0m\n");
     fprintf(stderr, "\033[32m\033[1mThis will be a depressingly long calculation ...     \033[0m");
   }

   #pragma omp parallel for private( arg, k, cosarg, sinarg)
   for (i=0;i<m;i++) {
      if ( m > 20000 )
        printf("%3d%%", (int)(100.0*done/m + 0.50) );
      arg = -(double)(6.2831853071795864769252867665590057683943387987) * (double)i / (double)m;
      for (k=0;k<m;k++) {
         cosarg = cos(k * arg);
         sinarg = sin(k * arg);
         x2[i] += (x1[k] * cosarg - y1[k] * sinarg);
         y2[i] += (x1[k] * sinarg + y1[k] * cosarg);
      }
      done++;
   }

   if ( m > 20000 )
     printf("");
     
   memcpy( x1, y1, m*sizeof(float) );
      
   for (i=0;i<m;i++) {
     y1[i] = 2*sqrt(x2[i]*x2[i]+y2[i]*y2[i])/m ;
   }

   data_DFT1 = x1;
   data_DFT2 = y1;
   
   free(x2);
   free(y2);
 
   HAVE_DFT = 2;

   return;
}



void myexit( code )
int code;
{
    if ( HAVE_TEMPFILE == YES )
        remove( tempfilename );

    exit( code );


}
