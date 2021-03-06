/*
 * main.c
 *
 *  Created on: 16.04.2014
 *      Author: jagiella
 */


#define ESC		27
#define SPACE	32
#define ENTER	10

#define SECONDS 1e6

#include <stdlib.h>

#include <ncurses.h>
#include <unistd.h>
#include <math.h>

#include "SDL2/SDL.h"
#include "SDL2_image/SDL_image.h"
#include "SDL2_ttf/SDL_ttf.h"

#include "tinyxml2.h"

void *
erealloc( void *v, size_t amt){
	if( v == 0){
		//fprintf(stderr, "malloc\n");
		v = malloc(amt);
	}else{
		//fprintf(stderr, "realloc\n");
		v = realloc(v, amt);
	}
    if( v == 0){
        fprintf(stderr, "out of mem\n");
        exit(EXIT_FAILURE);
    }
    return v;
}

char *mysprintf( char *filename, const char *fmt, int i){
	sprintf( filename, fmt, i);
	//fprintf( stderr, "Test %s\n", filename);
	return filename;
}

void flipHorizontally( SDL_Surface*& image );


class Spirit{
public:
	int            count_animations;
    int           *count_images;
    int            index_animation;
    int            index_images;
    SDL_Surface ***animations;

    SDL_Surface * spritesheet;
    SDL_Texture * texture;

    int            frames;
    SDL_Surface  **animation;

    SDL_Rect       rect;
    char           orientation;

    SDL_RendererFlip flip;

    enum { right, left, up, down};

    Spirit( SDL_Renderer * renderer){

    	// load animations


    	this->animations   = 0;
    	this->count_animations = 0;
    	this->count_images = 0;
    	this->index_images = 0;
    	this->index_animation = 0;

    	this->frames = 0;
    	this->animation = 0;


        //addAnimation( pattern);
    	/*spritesheet = IMG_Load("sprites-sheet-link.png");
    	SDL_SetColorKey( spritesheet, SDL_TRUE, SDL_MapRGB( spritesheet->format, 255, 255, 255 ) );
    	texture = SDL_CreateTextureFromSurface( renderer, spritesheet);*/

        // init position & size
    	this->rect.x=0;
    	this->rect.y=0;
    	this->rect.h=56;//44;
		this->rect.w=36;//38;
		//this->rect.h=48;//44;
		//this->rect.w=34;//38;
		this->orientation = right;
		flip = SDL_FLIP_NONE;

    };




    void addAnimation( const char *pattern){

    	int a = count_animations;
    	this->animations = (SDL_Surface***)erealloc( this->animations, sizeof(SDL_Surface**) * (a+1));
    	this->animations[a] = 0;

    	char filename[128];

		int i = 0;
		while( access( mysprintf( filename, pattern, i), F_OK ) != -1 )
		{
//			fprintf( stderr, "*\n");
			SDL_RWops *rwop	= SDL_RWFromFile( filename, "rb");
//	    	fprintf( stderr, "**\n");
	    	if( i == 0)
	    		this->animations[a] = (SDL_Surface**)erealloc( 0, sizeof(SDL_Surface*) * (i+10));
//	    	fprintf( stderr, "***\n");
			this->animations[a][i]	= IMG_LoadPNG_RW( rwop);
			if( this->animations[a]==0 )
				fprintf(stderr, "ERROR\n");
			i++;
		}

		count_animations ++;

		this->count_images = (int*)erealloc( this->count_images, sizeof(int) * (a+1));
		count_images[a] = i;

    };

    void loadFromFile( const char *pattern){

    	tinyxml2::XMLDocument doc;
    	doc.LoadFile( pattern );
    	doc.Print();

    	tinyxml2::XMLElement *animation = doc.FirstChildElement("Sprite")->FirstChildElement("Animation");
    	while( animation){
    		// NEW ANIMATION
    		int a = count_animations;
    		animations   = (SDL_Surface***)erealloc( this->animations, sizeof(SDL_Surface**) * (a+1));
    		count_images = (int*)erealloc( this->count_images, sizeof(int) * (a+1));
    		animations[a] = 0;
    		animations[a] = (SDL_Surface**)erealloc( 0, sizeof(SDL_Surface*) * (10));
			count_images[a] = 0;

    		tinyxml2::XMLElement *frame = animation->FirstChildElement("Frame");
    		int i=0;
    		while( frame){
    			// NEW FRAME
    			SDL_RWops *rwop	= SDL_RWFromFile( frame->FirstAttribute()->Value(), "rb");
    			animations[a][i] = IMG_LoadPNG_RW( rwop);
    			i ++;

    			frame = frame->NextSiblingElement();
    		};


    		count_images[a] = i;
    		count_animations ++;
    		animation = animation->NextSiblingElement();
    	};

    	//fprintf( stderr, (const char*) doc.FirstChildElement("Sprite")->FirstChildElement("Animation")->FirstChildElement("Frame")->FirstAttribute()->Value());


    	/*int a = count_animations;
    	this->animations = (SDL_Surface***)erealloc( this->animations, sizeof(SDL_Surface**) * (a+1));
    	this->animations[a] = 0;

    	char filename[128];

		int i = 0;
		while( access( mysprintf( filename, pattern, i), F_OK ) != -1 )
		{
//			fprintf( stderr, "*\n");
			SDL_RWops *rwop	= SDL_RWFromFile( filename, "rb");
//	    	fprintf( stderr, "**\n");
	    	if( i == 0)
	    		this->animations[a] = (SDL_Surface**)erealloc( 0, sizeof(SDL_Surface*) * (i+10));
//	    	fprintf( stderr, "***\n");
			this->animations[a][i]	= IMG_LoadPNG_RW( rwop);
			if( this->animations[a]==0 )
				fprintf(stderr, "ERROR\n");
			i++;
		}

		count_animations ++;

		this->count_images = (int*)erealloc( this->count_images, sizeof(int) * (a+1));
		count_images[a] = i;*/

    };

    void addFrames( const char *pattern){

        char filename[1024];

		int i = 0;
		while( access( mysprintf( filename, pattern, i), F_OK ) != -1 )
		{
	    	SDL_RWops *rwop	= SDL_RWFromFile( filename, "rb");
			this->animation = (SDL_Surface**)erealloc( this->animation, sizeof(SDL_Surface*) * (i+1));
			this->animation[i]	= IMG_LoadPNG_RW( rwop);
			i++;
		}

		this->frames = i;

		//exit(0);
    };

    void copyToRenderer( SDL_Renderer *renderer){
    	if( !count_animations) return;
    	//SDL_RWops   *rwops   = SDL_RWFromFile( "Zelda0.png", "rb");
    	//SDL_Surface *surface = IMG_LoadPNG_RW( rwops);
		 SDL_Texture *texture = SDL_CreateTextureFromSurface ( renderer, animations[index_animation][index_images]);
		 //SDL_RenderClear( renderer);
		 SDL_RenderCopyEx( renderer, texture, NULL, &rect, 0, 0, flip);
		 SDL_DestroyTexture(texture);

    };

 /*   void move( char new_orientation){
    	if( !count_animations) return;
    	index_images = (index_images+1)%count_images[index_animation]; // frames
    	//index_images = (index_images + 1) % count_images[index_animation];
    	if( orientation != new_orientation){
    		if( orientation == right || new_orientation == right) for( int i=0; i<count_images[0]; i++) flipHorizontally( animations[0][i] );

    		switch( new_orientation){
    		case left:	index_animation=0; flip = SDL_FLIP_HORIZONTAL; break;
    		case right: index_animation=0; flip = SDL_FLIP_NONE; break;
    		case up: 	index_animation=2; flip = SDL_FLIP_NONE; break;
    		case down:	index_animation=1; flip = SDL_FLIP_NONE; break;
    		}

        	orientation = new_orientation;
        	index_images = 0;
    	}
    	switch(orientation){
    	case right: rect.x+=5; break;
    	case left:  rect.x-=5; break;
    	case up:    rect.y-=5; break;
    	case down:  rect.y+=5; break;
    	}

    }*/

    void move( int h_axis, int v_axis){

    	if( h_axis || v_axis /*|| index_images*/ )
    		index_images = (index_images+1)%count_images[index_animation]; // frames
    	else
    		index_images = 0;

    	// change orientation?
    	char new_ori = orientation;
    	switch( h_axis ){
    	case -1: new_ori=left;  break;
    	case  0:  break;
    	case  1: new_ori=right; break;
    	}
    	switch( v_axis ){
    	case -1: new_ori=down; break;
    	case  0:  break;
    	case  1: new_ori=up; break;
    	}
    	if( orientation != new_ori){
    		//if( orientation == right || new_ori == right)


			switch( new_ori){
				case left:	index_animation=2; flip = SDL_FLIP_NONE; break;
				case right: index_animation=2; flip = SDL_FLIP_HORIZONTAL;  break;
				case up: 	index_animation=0; flip = SDL_FLIP_NONE;  break;
				case down:	index_animation=1; flip = SDL_FLIP_NONE;  break;
			}
    	}
    	orientation = new_ori;

    	int steps = sqrt( 50./(h_axis*h_axis + v_axis*v_axis) );
    	rect.x -= steps*h_axis;
    	rect.y -= steps*v_axis;

    }

};


void getSDLinput( int *dpad){

	SDL_Event event;

	if( SDL_PollEvent( &event ) )
			{
				//printf("key event");

				switch( event.type ){ /* Look for a keypress */
				case SDL_KEYDOWN:     /* Check the SDLKey values and move change the coords */
					switch( event.key.keysym.sym ){
					case SDLK_LEFT:		dpad[0] = 1; break;
					case SDLK_RIGHT:	dpad[1] = 1; break;
					case SDLK_UP:		dpad[2] = 1; break;
					case SDLK_DOWN:		dpad[3] = 1; break;
					case SDLK_ESCAPE:	exit(0); break;
					default: break;
					} break;
				case SDL_KEYUP:     /* Check the SDLKey values and move change the coords */
					switch( event.key.keysym.sym ){
					case SDLK_LEFT:		dpad[0] = 0; break;
					case SDLK_RIGHT:	dpad[1] = 0; break;
					case SDLK_UP:		dpad[2] = 0; break;
					case SDLK_DOWN:		dpad[3] = 0; break;
					default: break;
					} break;
				case SDL_JOYAXISMOTION:  /* Handle Joystick Motion */
					//if ( ( event.jaxis.value < -3200 ) || (event.jaxis.value > 3200 ) )
					{
						if( event.jaxis.axis == 0)
						{
							//mvprintw(2,2,"SDL_JOYAXISMOTION %+6i", event.jaxis.value);

							if( event.jaxis.value < -8000){
								dpad[0] = 1;
							}
							else if( event.jaxis.value > 8000)
								dpad[0] = -1;
							else
								dpad[0] = 0;

						            /* Left-right movement code goes here */
						}else

							if( event.jaxis.axis == 1)
							{
								//mvprintw(3,3,"SDL_JOYAXISMOTION %+6i", event.jaxis.value);

						            /* Up-Down movement code goes here */
								if( event.jaxis.value < -8000){
										dpad[2] = 1;
									}
									else if( event.jaxis.value > 8000)
										dpad[2] = -1;
									else
										dpad[2] = 0;

							}
					}//else{
						//dpad[0] = 1;
						//dpad[2] = 0;
					//}
				    break;
				}
			}
}

void updateJumpingBall( int *dpad, int &o, float *x, float *v, float *a, float *max)
{
	v[0] = v[0] - dpad[0] + dpad[1];
	v[1] = v[1] + (dpad[2] - dpad[3])*0.5;

	for( int d=0; d<2; d++){

		v[d] = 0.9*v[d] + 0.5*a[d];

		// REFLECTIVE BORDER
		if( x[d] < 1 || x[d] > max[d]){
			if( d==0)
				o = -o;
			v[d] = -v[d];
			x[d] = fmax(x[d],1);
			x[d] = fmin(x[d],max[d]);
		}


		// TORUS
		x[d] = x[d] + v[d];
		//x[d] = fmod( x[d] + v[d] + max[d], max[d]);

	}

}


void ncurses_printDpad( int *dpad, int offset_x, int offset_y){
	mvprintw(offset_y+0, offset_x, "+---+");
	mvprintw(offset_y+1, offset_x, "| %c |",                      (dpad[2]?'o':'.'));
	mvprintw(offset_y+2, offset_x, "|%c%c%c|", (dpad[0]?'o':'.'), (dpad[3]?'o':'.'), (dpad[1]?'o':'.'));
	mvprintw(offset_y+3, offset_x, "+---+");
}

void ncurses_printFrame( float *max){
    for( int ix=0;ix<max[0]; ix++){ // horizontal
   	 mvprintw(0, ix,   "-");
   	 mvprintw(max[1], ix,   "-");
    }
    for( int iy=0;iy<max[1]; iy++){ // vertical
   	 mvprintw(iy,0,    "|");
   	 mvprintw(iy,max[0],   "|");
    }

}

int main( int argc, char *argv[] )
{



	float max[2]={50, 40};
	float x[2] = {10,10};
	float v[2] = {0,0};
	float a[2] = {0,-1};
	int o = 1;

	int blcount=0;
	int blaster[100][3];

	int input=0;



    /* Initialise SDL */
     if( SDL_Init( SDL_INIT_VIDEO | SDL_INIT_JOYSTICK ) < 0){
         fprintf( stderr, "Could not initialise SDL: %s\n", SDL_GetError() );
         exit( -1 );
     }

     /* Set a video mode */
     SDL_Window *win = SDL_CreateWindow( "test",SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, 320, 200,
         		 SDL_WINDOW_INPUT_FOCUS);
     if( !win ){
         fprintf( stderr, "Could not set video mode: %s\n", SDL_GetError() );
         SDL_Quit();
         exit( -1 );
     }
     IMG_Init(IMG_INIT_PNG);

     TTF_Init();

     SDL_Renderer *ren = SDL_CreateRenderer(win, -1,
     	SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);

     // READ IMAGES
     Spirit *link = new Spirit( ren);

     link->loadFromFile( "LinkA.sprites" );
     //link->addFrames( "Zelda%d.png");
/*     link->addAnimation(  "ani_link_run_%d.png");
     link->addAnimation(  "ani_link_run_%d.png"); for( int i=0; i<link->count_images[0]; i++) flipHorizontally( link->animations[0][i] );
     link->addAnimation(  "ani_link_walk_front_%d.png");
     link->addAnimation(  "ani_link_walk_back_%d.png");*/


    // link->addFrames( "Zelda%d.png");
 /*         link->addAnimation(  "LinkA_walk_east%d.png");
          link->addAnimation(  "LinkA_walk_east%d.png"); for( int i=0; i<link->count_images[0]; i++) flipHorizontally( link->animations[0][i] );
          link->addAnimation(  "LinkA_walk_south%d.png");
          link->addAnimation(  "LinkA_walk_north%d.png");
*/
  		// Background
      	SDL_RWops *rwop	= SDL_RWFromFile( "dampes_house.png", "rb");
  		SDL_Surface *background	= IMG_LoadPNG_RW( rwop);
  		SDL_Texture *bg_texture = SDL_CreateTextureFromSurface ( ren, background);
  		SDL_Rect bg_rect;
  		bg_rect.x=0; bg_rect.y=0;bg_rect.w = 240; bg_rect.h = 160;

  		// Font
		int ptsize = 150;
		TTF_Font* myFont = TTF_OpenFontIndex( "returnofganon/ReturnofGanon.ttf", ptsize,0);
		SDL_Color fg; fg.r = 255; fg.b = 255; fg.g = 255;


     printf("%i joysticks were found.\n\n", SDL_NumJoysticks() );
     printf("The names of the joysticks are:\n");


/*     initscr();			// Start curses mode
     curs_set(0);
     nodelay(stdscr, TRUE);
     keypad(stdscr, TRUE);
     for( int i=0; i < SDL_NumJoysticks(); i++ )
     {
    	 mvprintw(1,1,"    %s\n", SDL_JoystickNameForIndex(i));
     }
     // paint border
     ncurses_printFrame( max);
*/

     SDL_JoystickEventState(SDL_ENABLE);
     SDL_Joystick *joystick = SDL_JoystickOpen(0);


	int dpad[4] = {0,0,0,0};

	int i=0; bool inc=true;
	while( 1 ){

		getSDLinput( dpad);

/*		ncurses_printDpad( dpad, 50, 0);
		mvprintw((int)(max[1]-x[1]), (int)x[0], " ");
		updateJumpingBall( dpad, o, x, v, a, max);
		mvprintw((int)(max[1]-x[1]), (int)x[0], "#");
		refresh();
*/





		// SDL ANIMATION
		/*if(dpad[1]==1)
			link->move( Spirit::right);
		if(dpad[0]==1)
			link->move( Spirit::left);
		if(dpad[2]==1)
			link->move( Spirit::up);
		if(dpad[3]==1)
			link->move( Spirit::down);
		 */

		// Clear old
		SDL_RenderClear( ren);

		// Background
		SDL_RenderCopy( ren, bg_texture, &bg_rect, 0);

		// Character
		link->move( dpad[0]-dpad[1], dpad[2]-dpad[3]);
		link->copyToRenderer( ren);

 /*       int width = 145/5;//32;
		SDL_Rect srcrect = { 20 + i*width, 80, width, 64 };
        SDL_Rect dstrect = { 10, 10, width, 64 };
 */       /*if( inc){
        	i++;
        	if(i==5){ i=4; inc=false;}
        }else{
        	i--;
        	if(i==0){ i=1; inc=true;}
        }*/
  //      i = (i+1)%6;


//	    SDL_RenderClear(  ren);
//	    SDL_RenderCopy(   ren, link->texture, &srcrect, &dstrect);

		// TEXT
		char dialogtext[512];
		if( link->rect.x > 100)
			sprintf(dialogtext, "Test!");
		else
			sprintf(dialogtext, "Down!");
		SDL_Surface *text = TTF_RenderText_Solid( myFont, dialogtext, fg);
		SDL_Texture *texture = SDL_CreateTextureFromSurface ( ren, text);
		int w, h; TTF_SizeText( myFont, dialogtext, &w, &h);
		SDL_Rect dstrect = { 10, 10, w, h };
		SDL_RenderCopy(   ren, texture, 0, &dstrect);
		//delete text;
		//delete texture;
		SDL_FreeSurface(text);
		SDL_DestroyTexture(texture);

		SDL_RenderPresent(ren);


		usleep( SECONDS/24. );
		 //usleep( 0.1*SECONDS );
	}


	//endwin();			/* End curses mode		  */

	return 0;
}


Uint32 get_pixel32( SDL_Surface *surface, int x, int y )
{
    //Convert the pixels to 32 bit
    Uint32 *pixels = (Uint32 *)surface->pixels;

    //Get the requested pixel
    return pixels[ ( y * surface->w ) + x ];
}



void put_pixel32( SDL_Surface *surface, int x, int y, Uint32 pixel )
{
    //Convert the pixels to 32 bit
    Uint32 *pixels = (Uint32 *)surface->pixels;

    //Set the pixel
    pixels[ ( y * surface->w ) + x ] = pixel;
}



void flipHorizontally( SDL_Surface*& image )
{
    // create a copy of the image
    SDL_Surface* flipped_image = SDL_CreateRGBSurface( SDL_SWSURFACE, image->w, image->h, image->format->BitsPerPixel,
        image->format->Rmask, image->format->Gmask, image->format->Bmask, image->format->Amask );

    // loop through pixels
    for( int y=0; y<image->h; y++ )
    {
        for( int x=0; x<image->w; x++ )
        {
            // copy pixels, but reverse the x pixels!
        	put_pixel32( flipped_image, x, y, get_pixel32(image, image->w - x - 1, y) );
        }
    }

    // free original and assign flipped to it
    SDL_FreeSurface( image );
    image = flipped_image;
}

