#ifndef GL_WINDOW_HPP_
#define GL_WINDOW_HPP_

#include <Agent.hpp>
#include <Boxes.hpp>
#include <Color.hpp>
#include <Macros.hpp>

#include <QtOpenGL/QGLWidget>
#include <QtOpenGL/QtOpenGL>

#include <GLUT/glut.h>

#include <QtGui/QMouseEvent>
#include <QtCore/QTime>


class WindowGL : public QGLWidget
{
    Q_OBJECT

public:
    WindowGL( AgentList<> *agentlist, Boxes<Agent*> *boxes, float *time, QWidget *parent = 0){ UNUSED(parent);
        QTime midnight(0, 0, 0);
        //qsrand(midnight.secsTo(QTime::currentTime()));

        this->time = time;
        image_idx = 0;
        al = agentlist;
        box = boxes;

        zoom = 0.;//-70;
        xRot = 0;
        yRot = 0;
        zRot = 0;

//        plotType = CellTypes;
        plotType = Pressure;
        //plotType = WNT;

        qtGreen = QColor::fromCmykF(0.40, 0.0, 1.0, 0.0);
        qtPurple = QColor::fromCmykF(0.39, 0.39, 0.0, 0.0);

        animationTimer.setSingleShot(false);
        connect(&animationTimer, SIGNAL(timeout()), this, SLOT(animate()));
        animationTimer.start(25);

        setAutoFillBackground(false);
        setMinimumSize(200, 200);
        this->resize(800,800);
        setWindowTitle(tr("Overpainting a Scene"));
        QIcon icon("animal_cell.jpg");
        setWindowIcon(icon);
    };
    ~WindowGL(){};

    void saveToFile ( QString saveFilename ){
    	QString saveExtension = "PNG";
    	int pos = saveFilename.lastIndexOf('.');
    	if (pos >= 0)
    	    saveExtension = saveFilename.mid(pos + 1);

    	this->grabFrameBuffer(this);
    	//glReadBuffer( GL_COLOR_BUFFER_BIT );

    	if(!QGLWidget::grabFrameBuffer( this).save(saveFilename, qPrintable(saveExtension)))
       	//if(!QPixmap::grabWidget( this).save(saveFilename, qPrintable(saveExtension)))
    	{
    	    // since you have a widget, just use grabWidget() here. winId() would possibly have
    	    // portability issues on other platforms.  qPrintable(saveExtension) is effectively
    	    // the same as saveExtension.toLocal8Bit().constData()

    	    //QMessageBox::warning(this, "File could not be saved", "ok", QMessageBox::Ok);
    	}
    	//QPixmap image = QPixmap::grabWidget( this );
    	//image.save( saveFilename, "PNG" );
    };

protected:
    void initializeGL(){
        // Set the clear color to black
        glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );

        // Set the drawing color to green
        glColor3f( 0.0f, 0.0f, 0.0f );

        // Light Source
        GLfloat light_ambient[] = { 1.0, 1.0, 1.0, 1.0 };
        GLfloat light_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
        GLfloat light_specular[] = { 1.0, 1.0, 1.0, 1.0 };
        GLfloat light_position[] = { box->X*box->latticeConstant/2, box->Y*box->latticeConstant/2, +100.0, 1.0 };

        glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
        glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
        glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
        glLightfv(GL_LIGHT0, GL_POSITION, light_position);

        // Light Options
        glEnable(GL_LIGHTING);
        glEnable(GL_LIGHT0);

        // Enable color tracking
        //glEnable( GL_COLOR_MATERIAL );
        // Set Material properties to follow glColor values
        glColorMaterial( GL_FRONT, GL_AMBIENT_AND_DIFFUSE );

        // Enable depth testing
        glEnable( GL_DEPTH_TEST );

        glShadeModel(GL_SMOOTH);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_BLEND);

    };
    void paintGL(){
    	makeCurrent();
    	//fprintf(stderr, "paint\n");
        // Clear the buffer with the current clearing color
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

        // Set the identity modelview transform
        glMatrixMode( GL_MODELVIEW );
        glLoadIdentity();


        // Translate the whole scene a little
        //glTranslatef( 0.0f, .0f, 0 );
//        glTranslatef( 0.0f, -75.0f, 0 );
        //glTranslatef( -100.0f, +100.0f, 200.0f );
        // CENTRALIZE
        //fprintf(stderr, "translate: %f %f\n", -box->X*box->latticeConstant/2., +box->Y*box->latticeConstant/2. );

//        glTranslatef( -box->X*box->latticeConstant/2., -box->Y*box->latticeConstant/2., 0.0f );
//        glRotatef( zoom, xRot, yRot, .0f );

        //glTranslatef( 0,100,0 );       // Translate away from the nucleus
        //glutSolidSphere(50, 15, 15 );

        //glFrustum(0, 200, 0, 200, 0, 200);

    	int delta = /*box->X*box->Y*box->Z*/+0;
    	int notch = /*box->X*box->Y*box->Z*/+1;


        if( plotType == Medium)
        for( int x=0; x<box->X; x++ )
        	for( int y=0; y<box->Y; y++ )
        		for( int z=0; z<box->Z; z++ ){
        			int dx=box->Y*box->Z;
        			int dy=box->Z;
        			int dz=1;
        			int j = x*dx + y*dy + z*dz;

        			glPushMatrix();
        			glTranslatef( (x+0.5)*box->latticeConstant,
        						  (y+0.5)*box->latticeConstant,
        						  (z-0.5)*box->latticeConstant );
        			//fprintf(stderr, "%e\n", box->molecules[j]);
        			GLfloat mat_diffuse[4] = { 1.f, 0.f, 0.f, box->molecules[j] };
        			glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
        			glutSolidCube( box->latticeConstant );
        			glPopMatrix();
        		}
        else

        for( int a=0; a<al->size(); a++)
        if(al->at(a)->position()[2] >= box->Z*box->latticeConstant/2.
        		|| (al->at(a)->type() != Agent::substrate &&  al->at(a)->_substrateContacts == 0))
        {
            glPushMatrix();                          // Save current modelview matrix

            //const float lightPos[] = { box->X*box->latticeConstant/2., box->Y*box->latticeConstant/2., +100.0f, 1.0f };
            //glLightfv( GL_LIGHT0, GL_POSITION, lightPos );

            float opacity = 1;
            //if(al->at(a)->position()[2] < 0)
            //	opacity = 0.1;

            GLfloat mat_diffuse[4] = { 0.f, 0.f, 0.f, opacity };

            switch( plotType){

            	case CellTypes:

                    switch( al->at(a)->type()){
                    case Agent::substrate: // blue
                    	mat_diffuse[2] = 0.8f;
                    	break;
                    case Agent::StemCell:  // green
                    	mat_diffuse[1] = 0.8f;
                    	break;
                    case Agent::Paneth:    // red
                    	mat_diffuse[0] = 0.8f;
                    	break;
                    case Agent::Goblet:    // light blue
                    	mat_diffuse[0] = 0.f;
                    	mat_diffuse[1] = 0.5f;
                    	mat_diffuse[2] = 0.8f;
                    	break;
                    case Agent::Enteroendocrine:    // light blue
                    	mat_diffuse[0] = 1.0f;
                    	mat_diffuse[1] = 1.0f;
                    	mat_diffuse[2] = 1.0f;
                    	break;
                    case Agent::Enterocyte:// yellow
                    	mat_diffuse[0] = 0.8f;
                    	mat_diffuse[1] = 0.8f;
                    	mat_diffuse[2] = 0.f;
                    	break;
                    case Agent::Absorptive:// orange
                    	mat_diffuse[0] = 1.0f;
                    	mat_diffuse[1] = 0.7f;
                    	mat_diffuse[2] = 0.f;
                    	break;
                    case Agent::Secretory:// Violet
                    	mat_diffuse[0] = 0.8f;
                    	mat_diffuse[1] = 0.f;
                    	mat_diffuse[2] = 0.8f;
                    	break;
                	}
            	break;

                case WNT:
                	mat_diffuse[0] = rgbformulae( 22, al->at(a)->wnt());
                	mat_diffuse[1] = rgbformulae( 13, al->at(a)->wnt());
                	mat_diffuse[2] = rgbformulae(-31, al->at(a)->wnt());

                	break;

                case Delta:
						mat_diffuse[0] = rgbformulae( 22, al->molecules[a*2+delta]);
						mat_diffuse[1] = rgbformulae( 13, al->molecules[a*2+delta]);
						mat_diffuse[2] = rgbformulae(-31, al->molecules[a*2+delta]);
                	break;
                case Notch:
                	if(al->at(a)->type() == Agent::substrate){
                		mat_diffuse[2] = 0.8f;
                	}else{
//						mat_diffuse[0] = rgbformulae( 22, al->at(a)->notch(al, box)/6.);
//						mat_diffuse[1] = rgbformulae( 13, al->at(a)->notch(al, box)/6.);
//						mat_diffuse[2] = rgbformulae(-31, al->at(a)->notch(al, box)/6.);
						mat_diffuse[0] = rgbformulae( 22, al->molecules[a*al->nbmolecules+plotMolecule]);
						mat_diffuse[1] = rgbformulae( 13, al->molecules[a*al->nbmolecules+plotMolecule]);
						mat_diffuse[2] = rgbformulae(-31, al->molecules[a*al->nbmolecules+plotMolecule]);
						//fprintf(stderr, "[ %e, %e ]\n", al->y[a*2]);
                	}
                	break;

                case Pressure:

                	//if(a==300)
                	//	fprintf(stdout, "\t\t\t\t\tpressure=%e (%p)\n", al->at(a)->pressure(), al->at(a));
                	mat_diffuse[0] = rgbformulae( 22, MINMAX( 0, fabs(al->at(a)->pressure())/al->pressureThreshold, 1));
                	mat_diffuse[1] = rgbformulae( 13, MINMAX( 0, fabs(al->at(a)->pressure())/al->pressureThreshold, 1));
                	mat_diffuse[2] = rgbformulae(-31, MINMAX( 0, fabs(al->at(a)->pressure())/al->pressureThreshold, 1));

                	break;

            }

//            if( al->at(a)->type() != Agent::substrate && al->at(a)->_substrateContacts == 0){
//				mat_diffuse[0] = 0;
//				mat_diffuse[1] = 0;
//				mat_diffuse[2] = 0;
//            }

            glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);

            //glRotatef( 10.0f, 0.0f, 1.0f, 0.0f );    // Make this orbit a little bit inclined to the basis vectors
            //glRotatef( m_theta3, 1.0f, 0.0f, 0.0f ); // Rotate by angle
            glTranslatef( al->at(a)->position()[0],
            		al->at(a)->position()[1],
            		al->at(a)->position()[2] );       // Translate away from the nucleus

            // LOOK TO CENTER
            /*glTranslatef( al->at(a)->position()[0]-box->X*box->latticeConstant/2,
            		-al->at(a)->position()[1]+box->X*box->latticeConstant/2,
            		al->at(a)->position()[2] );  */     // Translate away from the nucleus
            glutSolidSphere( al->at(a)->radius(), 15, 15 );         // Draw the electron

            glPopMatrix();                           // Restore the saved modelview matrix
        }

        //QPainter painter();
        //QFont font=painter.font() ;
         //font.setPointSize ( 18 );
         //font.setWeight(QFont::DemiBold);
         //QPainter::setFont(font);
        //QPainter painter(this);
        //painter.drawText(QPoint, "test");

        char ctext[512];
        // Plot Type
        switch(plotType){
        case CellTypes:
        	sprintf(ctext, "Color Code: Cell Types");
        	break;
        case Notch:
        	sprintf(ctext, "Color Code: Notch (%i)", plotMolecule);
    		break;
        case Delta:
        	sprintf(ctext, "Color Code: Delta");
    		break;
        case WNT:
        	sprintf(ctext, "Color Code: WNT");
        	break;
        case Pressure:
        	sprintf(ctext, "Color Code: Pressure");
        	break;
        }
    	renderText ( 50, 30, ctext, QFont("Arial", 30));

    	// Knockouts
        sprintf(ctext, "Knockout: ");
        if( Agent::WNTKnockout!=0)
        	sprintf(ctext, "%s WNT%c",ctext, (Agent::WNTKnockout==-1?'-':'+'));
        if( Agent::NotchKnockout)
        	sprintf(ctext, "%s Notch%c",ctext, (Agent::NotchKnockout==-1?'-':'+'));
        if( Agent::EphrinB12Knockout)
        	sprintf(ctext, "%s EphrinB1/B2%c",ctext, (Agent::EphrinB12Knockout==-1?'-':'+'));
        if( Agent::EphB2Knockout)
        	sprintf(ctext, "%s EphB2%c",ctext, (Agent::EphB2Knockout==-1?'-':'+'));
        if( Agent::EphB3Knockout)
        	sprintf(ctext, "%s EphB3%c",ctext, (Agent::EphB3Knockout==-1?'-':'+'));
        renderText ( 50, 60, ctext, QFont("Arial", 30));

        sprintf(ctext, "Time: %.2lf days", time[0]/24.);
        renderText ( 400, 30, ctext, QFont("Arial", 30));

    };
    void resizeGL(int w, int h){
        // Prevent a divde by zero
        if ( h == 0 )
            h = 1;

        // Set the viewport to window dimensions
        glViewport( 0, 0, w, h );

        // reset the coordinate system
        glMatrixMode( GL_PROJECTION );
        glLoadIdentity();

        // Establish the clipping volume by setting up an orthographic projection
        double aspectRatio = double( w ) / double( h );
        //double range = 300.0; //100.0;
        if ( true)//m_projectionMode == Orthogonal )
        {
            /*if ( w <=h )
                glOrtho( -range, range, -range / aspectRatio, range / aspectRatio, 2.0 * range, -2.0 * range );
            //glOrtho( 0, 2*range, -range / aspectRatio, range / aspectRatio, 2.0 * range, -2.0 * range );
            else
                glOrtho( -range * aspectRatio, range * aspectRatio, -range, range, 2.0 * range, -2.0 * range );
            //glOrtho( 0, 2*range * aspectRatio, -range, range, 2.0 * range, -2.0 * range );
*/
            /*glOrtho( 00, 200, // OK
            		-200, 00,
            		100, -100 );*/
            glOrtho( 0, box->X * box->latticeConstant, // OK
            		 0, box->Y * box->latticeConstant,
            		 0,-box->Z * box->latticeConstant );


        }
        else
        {
            double verticalViewingAngle = 60.0;
            gluPerspective( verticalViewingAngle, aspectRatio, 1.0, 200.0 );
        }

        glMatrixMode( GL_MODELVIEW );
        glLoadIdentity();
    };
    /*void mousePressEvent(QMouseEvent *event){
    	//event->
    };*/
    void mouseMoveEvent(QMouseEvent *event){
    	xRot = (float)event->pos().x();
    	yRot = (float)event->pos().y();
    	//fprintf(stderr, "rot = %e, %e\n", xRot,yRot);
    };

    void wheelEvent ( QWheelEvent * event ){
    	  //Add current step.
    	  //event->delta() can be negative or positive
    	zoom += event->delta()/120;

    	fprintf(stderr, "zoom = %e\n", zoom);
    	//this->setText("Total Steps: "+QString::number(zoom));
    };
    //void showEvent(QShowEvent *event){};
    void switchMod( char &x,char a, char b)
    {
    	//x = (x-a+1) % (b-a) + a;
    	x++;
    	if(x>b)
    		x=a;
    }
    void keyPressEvent(QKeyEvent * event){
    	switch( event->key()){
    	case Qt::Key_P:
    		plotType = Pressure;
    		break;

    	case Qt::Key_W:
    		plotType = WNT;
    		break;

    	case Qt::Key_T:
    		plotType = CellTypes;
    		break;

      	case Qt::Key_M:
       		plotType = Medium;
       		break;
      	case Qt::Key_N:
       		plotType = Notch;
       		plotMolecule = (plotMolecule + 1) % al->nbmolecules;
       		break;
      	case Qt::Key_D:
       		plotType = Delta;
       		break;
      	case Qt::Key_1:
      		switchMod( Agent::EphrinB12Knockout, -1, 1); break;
      	case Qt::Key_2:
      		switchMod( Agent::EphB2Knockout, -1, 1); break;
      	case Qt::Key_3:
      		switchMod( Agent::EphB3Knockout, -1, 1); break;
      	case Qt::Key_4:
       		switchMod( Agent::WNTKnockout, -1, 1); break;
      	case Qt::Key_5:
      		switchMod( Agent::NotchKnockout, -1, 1); break;
        	}
    };
private slots:
    void animate(){

        m_theta1 += 0.7f;
        if ( m_theta1 > 360.0f )
            m_theta1 -= 360.0f;

        m_theta2 += 1.5f;
        if ( m_theta2 > 360.0f )
            m_theta2 -= 360.0f;

        m_theta3 += 3.0f;
        if ( m_theta3 > 360.0f )
            m_theta3 -= 360.0f;

        updateGL();

        if(image_idx < (int)(*time*2.)){
      	  char filename[512];
      	  sprintf(filename, "test%03i.png", (int)(*time*2.));
      	  saveToFile(filename);
      	  image_idx++;
        }

   	    //update();
    };

private:
    //void createBubbles(int number){};
    //void drawInstructions(QPainter *painter){};
    //void setupViewport(int width, int height){};
//    ...
//    QtLogo *logo;
//    QList<Bubble*> bubbles;
    QTimer animationTimer;
    QColor qtGreen;
    QColor qtPurple;
	int xRot;
	int yRot;
	int zRot;
	float m_theta1,m_theta2,m_theta3;
	float zoom;

	AgentList<> *al;
	Boxes<Agent*> *box;

	float *time;
	int   image_idx;

	enum PlotTypes { CellTypes, WNT, Delta, Notch, Pressure, Medium};
	char plotType;
	int  plotMolecule;

};


#endif
