use tokio::net::TcpListener;
use tokio::io::AsyncReadExt;
use tokio::io::AsyncWriteExt;
use tokio::sync::Mutex;
use std::sync::{Arc};

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let listener = TcpListener::bind("127.0.0.1:55555").await?;
    let client = redis::Client::open("redis://db.occase.de/").unwrap();
    let con = Arc::new(Mutex::new(client.get_async_connection().await?));

    loop {
        let conn = Arc::clone(&con);
        let (mut socket, _) = listener.accept().await?;

        tokio::spawn(async move {
            let mut buf = [0; 1024];

            loop {
                let n = match socket.read(&mut buf).await {
                    Ok(n) if n == 0 => return,
                    Ok(n) => n,
                    Err(e) => {
                        eprintln!("failed to read from socket; err = {:?}", e);
                        return;
                    }
                };

                let mut local_conn = conn.lock().await;

                let result =
                    redis::cmd("PING")
                    .arg(&buf[0..n])
                    .query_async::<redis::aio::Connection, String>(&mut local_conn).await.unwrap();

                if let Err(e) = socket.write_all(result.as_bytes()).await {
                    eprintln!("failed to write to socket; err = {:?}", e);
                    return;
                }
            }
        });
    }
}
